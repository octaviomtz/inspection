import multiprocessing
import os
import pickle
import platform
import tarfile
import urllib.request
import warnings
from dataclasses import asdict, dataclass
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Literal, Optional

import click
import torch
from pytorch_lightning import Trainer, seed_everything
from pytorch_lightning.strategies import DDPStrategy
from pytorch_lightning.utilities import rank_zero_only
from rdkit import Chem
from tqdm import tqdm

from boltz.data import const
from boltz.data.module.inference import BoltzInferenceDataModule
from boltz.data.module.inferencev2 import Boltz2InferenceDataModule
from boltz.data.mol import load_canonicals
from boltz.data.msa.mmseqs2 import run_mmseqs2
from boltz.data.parse.a3m import parse_a3m
from boltz.data.parse.csv import parse_csv
from boltz.data.parse.fasta import parse_fasta
from boltz.data.parse.yaml import parse_yaml
from boltz.data.types import MSA, Manifest, Record
from boltz.data.write.writer import BoltzAffinityWriter, BoltzWriter
from boltz.model.models.boltz1 import Boltz1
from boltz.model.models.boltz2 import Boltz2

# Inspection hook (import at module level for optional use)
# The inspection module is at the repo root, not inside the package
# We need to add the repo root to sys.path to find it
import sys
_INSPECTION_AVAILABLE = False
_boltz_repo_root = Path(__file__).parent.parent.parent  # src/boltz/main.py -> repo root
_inspection_path = _boltz_repo_root / "inspection"
if _inspection_path.exists():
    sys.path.insert(0, str(_boltz_repo_root))
    try:
        from inspection.inspect_submodules import inspect_model_and_data
        from inspection.run_submodules import run_submodules_step_by_step
        from inspection.compare_outputs import compare_step_by_step_vs_forward
        from inspection.write_submodule_outputs import write_submodule_outputs
        _INSPECTION_AVAILABLE = True
    except ImportError:
        pass

CCD_URL = "https://huggingface.co/boltz-community/boltz-1/resolve/main/ccd.pkl"
MOL_URL = "https://huggingface.co/boltz-community/boltz-2/resolve/main/mols.tar"

BOLTZ1_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz1_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-1/resolve/main/boltz1_conf.ckpt",
]

BOLTZ2_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_conf.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_conf.ckpt",
]

BOLTZ2_AFFINITY_URL_WITH_FALLBACK = [
    "https://model-gateway.boltz.bio/boltz2_aff.ckpt",
    "https://huggingface.co/boltz-community/boltz-2/resolve/main/boltz2_aff.ckpt",
]


@dataclass
class BoltzProcessedInput:
    """Processed input data."""

    manifest: Manifest
    targets_dir: Path
    msa_dir: Path
    constraints_dir: Optional[Path] = None
    template_dir: Optional[Path] = None
    extra_mols_dir: Optional[Path] = None


@dataclass
class PairformerArgs:
    """Pairformer arguments."""

    num_blocks: int = 48
    num_heads: int = 16
    dropout: float = 0.0
    activation_checkpointing: bool = False
    offload_to_cpu: bool = False
    v2: bool = False


@dataclass
class PairformerArgsV2:
    """Pairformer arguments."""

    num_blocks: int = 64
    num_heads: int = 16
    dropout: float = 0.0
    activation_checkpointing: bool = False
    offload_to_cpu: bool = False
    v2: bool = True


@dataclass
class MSAModuleArgs:
    """MSA module arguments."""

    msa_s: int = 64
    msa_blocks: int = 4
    msa_dropout: float = 0.0
    z_dropout: float = 0.0
    use_paired_feature: bool = True
    pairwise_head_width: int = 32
    pairwise_num_heads: int = 4
    activation_checkpointing: bool = False
    offload_to_cpu: bool = False
    subsample_msa: bool = False
    num_subsampled_msa: int = 1024


@dataclass
class BoltzDiffusionParams:
    """Diffusion process parameters."""

    gamma_0: float = 0.605
    gamma_min: float = 1.107
    noise_scale: float = 0.901
    rho: float = 8
    step_scale: float = 1.638
    sigma_min: float = 0.0004
    sigma_max: float = 160.0
    sigma_data: float = 16.0
    P_mean: float = -1.2
    P_std: float = 1.5
    coordinate_augmentation: bool = True
    alignment_reverse_diff: bool = True
    synchronize_sigmas: bool = True
    use_inference_model_cache: bool = True


@dataclass
class Boltz2DiffusionParams:
    """Diffusion process parameters."""

    gamma_0: float = 0.8
    gamma_min: float = 1.0
    noise_scale: float = 1.003
    rho: float = 7
    step_scale: float = 1.5
    sigma_min: float = 0.0001
    sigma_max: float = 160.0
    sigma_data: float = 16.0
    P_mean: float = -1.2
    P_std: float = 1.5
    coordinate_augmentation: bool = True
    alignment_reverse_diff: bool = True
    synchronize_sigmas: bool = True


@dataclass
class BoltzSteeringParams:
    """Steering parameters."""

    fk_steering: bool = False
    num_particles: int = 3
    fk_lambda: float = 4.0
    fk_resampling_interval: int = 3
    physical_guidance_update: bool = False
    contact_guidance_update: bool = True
    cdr3_steering: bool = False
    antigen_steering: bool = False
    num_gd_steps: int = 20


@rank_zero_only
def download_boltz1(cache: Path) -> None:
    """Download all the required data.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    # Download CCD
    ccd = cache / "ccd.pkl"
    if not ccd.exists():
        click.echo(
            f"Downloading the CCD dictionary to {ccd}. You may "
            "change the cache directory with the --cache flag."
        )
        urllib.request.urlretrieve(CCD_URL, str(ccd))  # noqa: S310

    # Download model
    model = cache / "boltz1_conf.ckpt"
    if not model.exists():
        click.echo(
            f"Downloading the model weights to {model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ1_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(model))  # noqa: S310
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ1_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue


@rank_zero_only
def download_boltz2(cache: Path) -> None:
    """Download all the required data.

    Parameters
    ----------
    cache : Path
        The cache directory.

    """
    # Download CCD
    mols = cache / "mols"
    tar_mols = cache / "mols.tar"
    if not tar_mols.exists():
        click.echo(
            f"Downloading the CCD data to {tar_mols}. "
            "This may take a bit of time. You may change the cache directory "
            "with the --cache flag."
        )
        urllib.request.urlretrieve(MOL_URL, str(tar_mols))  # noqa: S310
    if not mols.exists():
        click.echo(
            f"Extracting the CCD data to {mols}. "
            "This may take a bit of time. You may change the cache directory "
            "with the --cache flag."
        )
        with tarfile.open(str(tar_mols), "r") as tar:
            tar.extractall(cache)  # noqa: S202

    # Download model
    model = cache / "boltz2_conf.ckpt"
    if not model.exists():
        click.echo(
            f"Downloading the Boltz-2 weights to {model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ2_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(model))  # noqa: S310
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue

    # Download affinity model
    affinity_model = cache / "boltz2_aff.ckpt"
    if not affinity_model.exists():
        click.echo(
            f"Downloading the Boltz-2 affinity weights to {affinity_model}. You may "
            "change the cache directory with the --cache flag."
        )
        for i, url in enumerate(BOLTZ2_AFFINITY_URL_WITH_FALLBACK):
            try:
                urllib.request.urlretrieve(url, str(affinity_model))  # noqa: S310
                break
            except Exception as e:  # noqa: BLE001
                if i == len(BOLTZ2_AFFINITY_URL_WITH_FALLBACK) - 1:
                    msg = f"Failed to download model from all URLs. Last error: {e}"
                    raise RuntimeError(msg) from e
                continue


def get_cache_path() -> str:
    """Determine the cache path, prioritising the BOLTZ_CACHE environment variable.

    Returns
    -------
    str: Path
        Path to use for boltz cache location.

    """
    env_cache = os.environ.get("BOLTZ_CACHE")
    if env_cache:
        resolved_cache = Path(env_cache).expanduser().resolve()
        if not resolved_cache.is_absolute():
            msg = f"BOLTZ_CACHE must be an absolute path, got: {env_cache}"
            raise ValueError(msg)
        return str(resolved_cache)

    return str(Path("~/.boltz").expanduser())


def check_inputs(data: Path) -> list[Path]:
    """Check the input data and output directory.

    Parameters
    ----------
    data : Path
        The input data.

    Returns
    -------
    list[Path]
        The list of input data.

    """
    click.echo("Checking input data.")

    # Check if data is a directory
    if data.is_dir():
        data: list[Path] = list(data.glob("*"))

        # Filter out non .fasta or .yaml files, raise
        # an error on directory and other file types
        for d in data:
            if d.is_dir():
                msg = f"Found directory {d} instead of .fasta or .yaml."
                raise RuntimeError(msg)
            if d.suffix.lower() not in (".fa", ".fas", ".fasta", ".yml", ".yaml"):
                msg = (
                    f"Unable to parse filetype {d.suffix}, "
                    "please provide a .fasta or .yaml file."
                )
                raise RuntimeError(msg)
    else:
        data = [data]

    return data


def filter_inputs_structure(
    manifest: Manifest,
    outdir: Path,
    override: bool = False,
) -> Manifest:
    """Filter the manifest to only include missing predictions.

    Parameters
    ----------
    manifest : Manifest
        The manifest of the input data.
    outdir : Path
        The output directory.
    override: bool
        Whether to override existing predictions.

    Returns
    -------
    Manifest
        The manifest of the filtered input data.

    """
    # Check if existing predictions are found (only top-level prediction folders)
    pred_dir = outdir / "predictions"
    if pred_dir.exists():
        existing = {d.name for d in pred_dir.iterdir() if d.is_dir()}
    else:
        existing = set()

    # Remove them from the input data
    if existing and not override:
        manifest = Manifest([r for r in manifest.records if r.id not in existing])
        msg = (
            f"Found some existing predictions ({len(existing)}), "
            f"skipping and running only the missing ones, "
            "if any. If you wish to override these existing "
            "predictions, please set the --override flag."
        )
        click.echo(msg)
    elif existing and override:
        msg = f"Found {len(existing)} existing predictions, will override."
        click.echo(msg)

    return manifest


def filter_inputs_affinity(
    manifest: Manifest,
    outdir: Path,
    override: bool = False,
) -> Manifest:
    """Check the input data and output directory for affinity.

    Parameters
    ----------
    manifest : Manifest
        The manifest.
    outdir : Path
        The output directory.
    override: bool
        Whether to override existing predictions.

    Returns
    -------
    Manifest
        The manifest of the filtered input data.

    """
    click.echo("Checking input data for affinity.")

    # Get all affinity targets
    existing = {
        r.id
        for r in manifest.records
        if r.affinity
        and (outdir / "predictions" / r.id / f"affinity_{r.id}.json").exists()
    }

    # Remove them from the input data
    if existing and not override:
        manifest = Manifest([r for r in manifest.records if r.id not in existing])
        num_skipped = len(existing)
        msg = (
            f"Found some existing affinity predictions ({num_skipped}), "
            f"skipping and running only the missing ones, "
            "if any. If you wish to override these existing "
            "affinity predictions, please set the --override flag."
        )
        click.echo(msg)
    elif existing and override:
        msg = "Found existing affinity predictions, will override."
        click.echo(msg)

    return manifest


def compute_msa(
    data: dict[str, str],
    target_id: str,
    msa_dir: Path,
    msa_server_url: str,
    msa_pairing_strategy: str,
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    api_key_header: Optional[str] = None,
    api_key_value: Optional[str] = None,
) -> None:
    """Compute the MSA for the input data.

    Parameters
    ----------
    data : dict[str, str]
        The input protein sequences.
    target_id : str
        The target id.
    msa_dir : Path
        The msa directory.
    msa_server_url : str
        The MSA server URL.
    msa_pairing_strategy : str
        The MSA pairing strategy.
    msa_server_username : str, optional
        Username for basic authentication with MSA server.
    msa_server_password : str, optional
        Password for basic authentication with MSA server.
    api_key_header : str, optional
        Custom header key for API key authentication (default: X-API-Key).
    api_key_value : str, optional
        Custom header value for API key authentication (overrides --api_key if set).

    """
    click.echo(f"Calling MSA server for target {target_id} with {len(data)} sequences")
    click.echo(f"MSA server URL: {msa_server_url}")
    click.echo(f"MSA pairing strategy: {msa_pairing_strategy}")
    
    # Construct auth headers if API key header/value is provided
    auth_headers = None
    if api_key_value:
        key = api_key_header if api_key_header else "X-API-Key"
        value = api_key_value
        auth_headers = {
            "Content-Type": "application/json",
            key: value
        }
        click.echo(f"Using API key authentication for MSA server (header: {key})")
    elif msa_server_username and msa_server_password:
        click.echo("Using basic authentication for MSA server")
    else:
        click.echo("No authentication provided for MSA server")
    
    if len(data) > 1:
        paired_msas = run_mmseqs2(
            list(data.values()),
            msa_dir / f"{target_id}_paired_tmp",
            use_env=True,
            use_pairing=True,
            host_url=msa_server_url,
            pairing_strategy=msa_pairing_strategy,
            msa_server_username=msa_server_username,
            msa_server_password=msa_server_password,
            auth_headers=auth_headers,
        )
    else:
        paired_msas = [""] * len(data)

    unpaired_msa = run_mmseqs2(
        list(data.values()),
        msa_dir / f"{target_id}_unpaired_tmp",
        use_env=True,
        use_pairing=False,
        host_url=msa_server_url,
        pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        auth_headers=auth_headers,
    )

    for idx, name in enumerate(data):
        # Get paired sequences
        paired = paired_msas[idx].strip().splitlines()
        paired = paired[1::2]  # ignore headers
        paired = paired[: const.max_paired_seqs]

        # Set key per row and remove empty sequences
        keys = [idx for idx, s in enumerate(paired) if s != "-" * len(s)]
        paired = [s for s in paired if s != "-" * len(s)]

        # Combine paired-unpaired sequences
        unpaired = unpaired_msa[idx].strip().splitlines()
        unpaired = unpaired[1::2]
        unpaired = unpaired[: (const.max_msa_seqs - len(paired))]
        if paired:
            unpaired = unpaired[1:]  # ignore query is already present

        # Combine
        seqs = paired + unpaired
        keys = keys + [-1] * len(unpaired)

        # Dump MSA
        csv_str = ["key,sequence"] + [f"{key},{seq}" for key, seq in zip(keys, seqs)]

        msa_path = msa_dir / f"{name}.csv"
        with msa_path.open("w") as f:
            f.write("\n".join(csv_str))


def process_input(  # noqa: C901, PLR0912, PLR0915, D103
    path: Path,
    ccd: dict,
    msa_dir: Path,
    mol_dir: Path,
    boltz2: bool,
    use_msa_server: bool,
    msa_server_url: str,
    msa_pairing_strategy: str,
    msa_server_username: Optional[str],
    msa_server_password: Optional[str],
    api_key_header: Optional[str],
    api_key_value: Optional[str],
    max_msa_seqs: int,
    processed_msa_dir: Path,
    processed_constraints_dir: Path,
    processed_templates_dir: Path,
    processed_mols_dir: Path,
    structure_dir: Path,
    records_dir: Path,
) -> None:
    try:
        # Parse data
        if path.suffix.lower() in (".fa", ".fas", ".fasta"):
            target = parse_fasta(path, ccd, mol_dir, boltz2)
        elif path.suffix.lower() in (".yml", ".yaml"):
            target = parse_yaml(path, ccd, mol_dir, boltz2)
        elif path.is_dir():
            msg = f"Found directory {path} instead of .fasta or .yaml, skipping."
            raise RuntimeError(msg)  # noqa: TRY301
        else:
            msg = (
                f"Unable to parse filetype {path.suffix}, "
                "please provide a .fasta or .yaml file."
            )
            raise RuntimeError(msg)  # noqa: TRY301

        # Get target id
        target_id = target.record.id

        # Get all MSA ids and decide whether to generate MSA
        to_generate = {}
        prot_id = const.chain_type_ids["PROTEIN"]
        for chain in target.record.chains:
            # Add to generate list, assigning entity id
            if (chain.mol_type == prot_id) and (chain.msa_id == 0):
                entity_id = chain.entity_id
                msa_id = f"{target_id}_{entity_id}"
                to_generate[msa_id] = target.sequences[entity_id]
                chain.msa_id = msa_dir / f"{msa_id}.csv"

            # We do not support msa generation for non-protein chains
            elif chain.msa_id == 0:
                chain.msa_id = -1

        # Generate MSA
        if to_generate and not use_msa_server:
            msg = "Missing MSA's in input and --use_msa_server flag not set."
            raise RuntimeError(msg)  # noqa: TRY301

        if to_generate:
            msg = f"Generating MSA for {path} with {len(to_generate)} protein entities."
            click.echo(msg)
            compute_msa(
                data=to_generate,
                target_id=target_id,
                msa_dir=msa_dir,
                msa_server_url=msa_server_url,
                msa_pairing_strategy=msa_pairing_strategy,
                msa_server_username=msa_server_username,
                msa_server_password=msa_server_password,
                api_key_header=api_key_header,
                api_key_value=api_key_value,
            )

        # Parse MSA data
        msas = sorted({c.msa_id for c in target.record.chains if c.msa_id != -1})
        msa_id_map = {}
        for msa_idx, msa_id in enumerate(msas):
            # Check that raw MSA exists
            msa_path = Path(msa_id)
            if not msa_path.exists():
                msg = f"MSA file {msa_path} not found."
                raise FileNotFoundError(msg)  # noqa: TRY301

            # Dump processed MSA
            processed = processed_msa_dir / f"{target_id}_{msa_idx}.npz"
            msa_id_map[msa_id] = f"{target_id}_{msa_idx}"
            if not processed.exists():
                # Parse A3M
                if msa_path.suffix == ".a3m":
                    msa: MSA = parse_a3m(
                        msa_path,
                        taxonomy=None,
                        max_seqs=max_msa_seqs,
                    )
                elif msa_path.suffix == ".csv":
                    msa: MSA = parse_csv(msa_path, max_seqs=max_msa_seqs)
                else:
                    msg = f"MSA file {msa_path} not supported, only a3m or csv."
                    raise RuntimeError(msg)  # noqa: TRY301

                msa.dump(processed)

        # Modify records to point to processed MSA
        for c in target.record.chains:
            if (c.msa_id != -1) and (c.msa_id in msa_id_map):
                c.msa_id = msa_id_map[c.msa_id]

        # Dump templates
        for template_id, template in target.templates.items():
            name = f"{target.record.id}_{template_id}.npz"
            template_path = processed_templates_dir / name
            template.dump(template_path)

        # Dump constraints
        constraints_path = processed_constraints_dir / f"{target.record.id}.npz"
        target.residue_constraints.dump(constraints_path)

        # Dump extra molecules
        Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)
        with (processed_mols_dir / f"{target.record.id}.pkl").open("wb") as f:
            pickle.dump(target.extra_mols, f)

        # Dump structure
        struct_path = structure_dir / f"{target.record.id}.npz"
        target.structure.dump(struct_path)

        # Dump record
        record_path = records_dir / f"{target.record.id}.json"
        target.record.dump(record_path)

    except Exception as e:  # noqa: BLE001
        import traceback

        traceback.print_exc()
        print(f"Failed to process {path}. Skipping. Error: {e}.")  # noqa: T201


@rank_zero_only
def process_inputs(
    data: list[Path],
    out_dir: Path,
    ccd_path: Path,
    mol_dir: Path,
    msa_server_url: str,
    msa_pairing_strategy: str,
    max_msa_seqs: int = 8192,
    use_msa_server: bool = False,
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    api_key_header: Optional[str] = None,
    api_key_value: Optional[str] = None,
    boltz2: bool = False,
    preprocessing_threads: int = 1,
) -> Manifest:
    """Process the input data and output directory.

    Parameters
    ----------
    data : list[Path]
        The input data.
    out_dir : Path
        The output directory.
    ccd_path : Path
        The path to the CCD dictionary.
    max_msa_seqs : int, optional
        Max number of MSA sequences, by default 8192.
    use_msa_server : bool, optional
        Whether to use the MMSeqs2 server for MSA generation, by default False.
    msa_server_username : str, optional
        Username for basic authentication with MSA server, by default None.
    msa_server_password : str, optional
        Password for basic authentication with MSA server, by default None.
    api_key_header : str, optional
        Custom header key for API key authentication (default: X-API-Key).
    api_key_value : str, optional
        Custom header value for API key authentication (overrides --api_key if set).
    boltz2: bool, optional
        Whether to use Boltz2, by default False.
    preprocessing_threads: int, optional
        The number of threads to use for preprocessing, by default 1.

    Returns
    -------
    Manifest
        The manifest of the processed input data.

    """
    # Validate mutually exclusive authentication methods
    has_basic_auth = msa_server_username and msa_server_password
    has_api_key = api_key_value is not None
    
    if has_basic_auth and has_api_key:
        raise ValueError(
            "Cannot use both basic authentication (--msa_server_username/--msa_server_password) "
            "and API key authentication (--api_key_header/--api_key_value). Please use only one authentication method."
        )

    # Check if records exist at output path
    records_dir = out_dir / "processed" / "records"
    if records_dir.exists():
        # Load existing records
        existing = [Record.load(p) for p in records_dir.glob("*.json")]
        processed_ids = {record.id for record in existing}

        # Filter to missing only
        data = [d for d in data if d.stem not in processed_ids]

        # Nothing to do, update the manifest and return
        if data:
            click.echo(
                f"Found {len(existing)} existing processed inputs, skipping them."
            )
        else:
            click.echo("All inputs are already processed.")
            updated_manifest = Manifest(existing)
            updated_manifest.dump(out_dir / "processed" / "manifest.json")

    # Create output directories
    msa_dir = out_dir / "msa"
    records_dir = out_dir / "processed" / "records"
    structure_dir = out_dir / "processed" / "structures"
    processed_msa_dir = out_dir / "processed" / "msa"
    processed_constraints_dir = out_dir / "processed" / "constraints"
    processed_templates_dir = out_dir / "processed" / "templates"
    processed_mols_dir = out_dir / "processed" / "mols"
    predictions_dir = out_dir / "predictions"

    out_dir.mkdir(parents=True, exist_ok=True)
    msa_dir.mkdir(parents=True, exist_ok=True)
    records_dir.mkdir(parents=True, exist_ok=True)
    structure_dir.mkdir(parents=True, exist_ok=True)
    processed_msa_dir.mkdir(parents=True, exist_ok=True)
    processed_constraints_dir.mkdir(parents=True, exist_ok=True)
    processed_templates_dir.mkdir(parents=True, exist_ok=True)
    processed_mols_dir.mkdir(parents=True, exist_ok=True)
    predictions_dir.mkdir(parents=True, exist_ok=True)

    # Load CCD
    if boltz2:
        ccd = load_canonicals(mol_dir)
    else:
        with ccd_path.open("rb") as file:
            ccd = pickle.load(file)  # noqa: S301

    # Create partial function
    process_input_partial = partial(
        process_input,
        ccd=ccd,
        msa_dir=msa_dir,
        mol_dir=mol_dir,
        boltz2=boltz2,
        use_msa_server=use_msa_server,
        msa_server_url=msa_server_url,
        msa_pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        api_key_header=api_key_header,
        api_key_value=api_key_value,
        max_msa_seqs=max_msa_seqs,
        processed_msa_dir=processed_msa_dir,
        processed_constraints_dir=processed_constraints_dir,
        processed_templates_dir=processed_templates_dir,
        processed_mols_dir=processed_mols_dir,
        structure_dir=structure_dir,
        records_dir=records_dir,
    )

    # Parse input data
    preprocessing_threads = min(preprocessing_threads, len(data))
    click.echo(f"Processing {len(data)} inputs with {preprocessing_threads} threads.")

    if preprocessing_threads > 1 and len(data) > 1:
        with Pool(preprocessing_threads) as pool:
            list(tqdm(pool.imap(process_input_partial, data), total=len(data)))
    else:
        for path in tqdm(data):
            process_input_partial(path)

    # Load all records and write manifest
    records = [Record.load(p) for p in records_dir.glob("*.json")]
    manifest = Manifest(records)
    manifest.dump(out_dir / "processed" / "manifest.json")


@click.group()
def cli() -> None:
    """Boltz."""
    return


@cli.command()
@click.argument("data", type=click.Path(exists=True))
@click.option(
    "--out_dir",
    type=click.Path(exists=False),
    help="The path where to save the predictions.",
    default="./",
)
@click.option(
    "--cache",
    type=click.Path(exists=False),
    help=(
        "The directory where to download the data and model. "
        "Default is ~/.boltz, or $BOLTZ_CACHE if set."
    ),
    default=get_cache_path,
)
@click.option(
    "--checkpoint",
    type=click.Path(exists=True),
    help="An optional checkpoint, will use the provided Boltz-1 model by default.",
    default=None,
)
@click.option(
    "--devices",
    type=int,
    help="The number of devices to use for prediction. Default is 1.",
    default=1,
)
@click.option(
    "--accelerator",
    type=click.Choice(["gpu", "cpu", "tpu"]),
    help="The accelerator to use for prediction. Default is gpu.",
    default="gpu",
)
@click.option(
    "--recycling_steps",
    type=int,
    help="The number of recycling steps to use for prediction. Default is 3.",
    default=3,
)
@click.option(
    "--sampling_steps",
    type=int,
    help="The number of sampling steps to use for prediction. Default is 200.",
    default=200,
)
@click.option(
    "--diffusion_samples",
    type=int,
    help="The number of diffusion samples to use for prediction. Default is 1.",
    default=1,
)
@click.option(
    "--max_parallel_samples",
    type=int,
    help="The maximum number of samples to predict in parallel. Default is None.",
    default=5,
)
@click.option(
    "--step_scale",
    type=float,
    help=(
        "The step size is related to the temperature at "
        "which the diffusion process samples the distribution. "
        "The lower the higher the diversity among samples "
        "(recommended between 1 and 2). "
        "Default is 1.638 for Boltz-1 and 1.5 for Boltz-2. "
        "If not provided, the default step size will be used."
    ),
    default=None,
)
@click.option(
    "--write_full_pae",
    type=bool,
    is_flag=True,
    help="Whether to dump the pae into a npz file. Default is True.",
)
@click.option(
    "--write_full_pde",
    type=bool,
    is_flag=True,
    help="Whether to dump the pde into a npz file. Default is False.",
)
@click.option(
    "--output_format",
    type=click.Choice(["pdb", "mmcif"]),
    help="The output format to use for the predictions. Default is mmcif.",
    default="mmcif",
)
@click.option(
    "--num_workers",
    type=int,
    help="The number of dataloader workers to use for prediction. Default is 2.",
    default=2,
)
@click.option(
    "--override",
    is_flag=True,
    help="Whether to override existing found predictions. Default is False.",
)
@click.option(
    "--seed",
    type=int,
    help="Seed to use for random number generator. Default is None (no seeding).",
    default=None,
)
@click.option(
    "--use_msa_server",
    is_flag=True,
    help="Whether to use the MMSeqs2 server for MSA generation. Default is False.",
)
@click.option(
    "--msa_server_url",
    type=str,
    help="MSA server url. Used only if --use_msa_server is set. ",
    default="https://api.colabfold.com",
)
@click.option(
    "--msa_pairing_strategy",
    type=str,
    help=(
        "Pairing strategy to use. Used only if --use_msa_server is set. "
        "Options are 'greedy' and 'complete'"
    ),
    default="greedy",
)
@click.option(
    "--msa_server_username",
    type=str,
    help="MSA server username for basic auth. Used only if --use_msa_server is set. Can also be set via BOLTZ_MSA_USERNAME environment variable.",
    default=None,
)
@click.option(
    "--msa_server_password",
    type=str,
    help="MSA server password for basic auth. Used only if --use_msa_server is set. Can also be set via BOLTZ_MSA_PASSWORD environment variable.",
    default=None,
)
@click.option(
    "--api_key_header",
    type=str,
    help="Custom header key for API key authentication (default: X-API-Key).",
    default=None,
)
@click.option(
    "--api_key_value",
    type=str,
    help="Custom header value for API key authentication.",
    default=None,
)
@click.option(
    "--use_potentials",
    is_flag=True,
    help="Whether to use potentials for steering. Default is False.",
)
@click.option(
    "--cdr3_steering",
    is_flag=True,
    help="Enable CDR3 conformation steering via backbone dihedral angle potentials. "
         "Requires --use_potentials. Define CDR3 regions using cdr3_conformation constraints in YAML.",
)
@click.option(
    "--antigen_steering",
    is_flag=True,
    help="Enable antigen orientation steering to maximize CDR-antigen contact. "
         "Requires --use_potentials. Define using antigen_orientation constraints in YAML.",
)
@click.option(
    "--num_particles",
    type=int,
    default=3,
    help="Number of particles for Feynman-Kac steering. Lower values reduce memory usage. "
         "With steering enabled, total samples = diffusion_samples × num_particles. Default is 3.",
)
@click.option(
    "--model",
    default="boltz2",
    type=click.Choice(["boltz1", "boltz2"]),
    help="The model to use for prediction. Default is boltz2.",
)
@click.option(
    "--method",
    type=str,
    help="The method to use for prediction. Default is None.",
    default=None,
)
@click.option(
    "--preprocessing-threads",
    type=int,
    help="The number of threads to use for preprocessing. Default is 1.",
    default=multiprocessing.cpu_count(),
)
@click.option(
    "--affinity_mw_correction",
    is_flag=True,
    type=bool,
    help="Whether to add the Molecular Weight correction to the affinity value head.",
)
@click.option(
    "--sampling_steps_affinity",
    type=int,
    help="The number of sampling steps to use for affinity prediction. Default is 200.",
    default=200,
)
@click.option(
    "--diffusion_samples_affinity",
    type=int,
    help="The number of diffusion samples to use for affinity prediction. Default is 5.",
    default=5,
)
@click.option(
    "--affinity_checkpoint",
    type=click.Path(exists=True),
    help="An optional checkpoint, will use the provided Boltz-1 model by default.",
    default=None,
)
@click.option(
    "--max_msa_seqs",
    type=int,
    help="The maximum number of MSA sequences to use for prediction. Default is 8192.",
    default=8192,
)
@click.option(
    "--subsample_msa",
    is_flag=True,
    help="Whether to subsample the MSA. Default is True.",
)
@click.option(
    "--num_subsampled_msa",
    type=int,
    help="The number of MSA sequences to subsample. Default is 1024.",
    default=1024,
)
@click.option(
    "--no_kernels",
    is_flag=True,
    help="Whether to disable the kernels. Default False",
)
@click.option(
    "--write_embeddings",
    is_flag=True,
    help=" to dump the s and z embeddings into a npz file. Default is False.",
)
@click.option(
    "--step_by_step",
    is_flag=True,
    help="Run the model step by step, saving intermediate submodule outputs. Default is False.",
)
@click.option(
    "--canonical_ensemble",
    is_flag=True,
    help="Run canonical ensemble sweep over beta values. Requires canonical_ensemble constraint in YAML.",
)
def predict(  # noqa: C901, PLR0915, PLR0912
    data: str,
    out_dir: str,
    cache: str = "~/.boltz",
    checkpoint: Optional[str] = None,
    affinity_checkpoint: Optional[str] = None,
    devices: int = 1,
    accelerator: str = "gpu",
    recycling_steps: int = 3,
    sampling_steps: int = 200,
    diffusion_samples: int = 1,
    sampling_steps_affinity: int = 200,
    diffusion_samples_affinity: int = 3,
    max_parallel_samples: Optional[int] = None,
    step_scale: Optional[float] = None,
    write_full_pae: bool = False,
    write_full_pde: bool = False,
    output_format: Literal["pdb", "mmcif"] = "mmcif",
    num_workers: int = 2,
    override: bool = False,
    seed: Optional[int] = None,
    use_msa_server: bool = False,
    msa_server_url: str = "https://api.colabfold.com",
    msa_pairing_strategy: str = "greedy",
    msa_server_username: Optional[str] = None,
    msa_server_password: Optional[str] = None,
    api_key_header: Optional[str] = None,
    api_key_value: Optional[str] = None,
    use_potentials: bool = False,
    cdr3_steering: bool = False,
    antigen_steering: bool = False,
    num_particles: int = 3,
    model: Literal["boltz1", "boltz2"] = "boltz2",
    method: Optional[str] = None,
    affinity_mw_correction: Optional[bool] = False,
    preprocessing_threads: int = 1,
    max_msa_seqs: int = 8192,
    subsample_msa: bool = True,
    num_subsampled_msa: int = 1024,
    no_kernels: bool = False,
    write_embeddings: bool = False,
    step_by_step: bool = False,
    canonical_ensemble: bool = False,
) -> None:
    """Run predictions with Boltz."""
    # If cpu, write a friendly warning
    if accelerator == "cpu":
        msg = "Running on CPU, this will be slow. Consider using a GPU."
        click.echo(msg)

    # Supress some lightning warnings
    warnings.filterwarnings(
        "ignore", ".*that has Tensor Cores. To properly utilize them.*"
    )

    # Set no grad
    torch.set_grad_enabled(False)

    # Ignore matmul precision warning
    torch.set_float32_matmul_precision("highest")

    # Set rdkit pickle logic
    Chem.SetDefaultPickleProperties(Chem.PropertyPickleOptions.AllProps)

    # Set seed if desired
    if seed is not None:
        seed_everything(seed)

    for key in ["CUEQ_DEFAULT_CONFIG", "CUEQ_DISABLE_AOT_TUNING"]:
        # Disable kernel tuning by default,
        # but do not modify envvar if already set by caller
        os.environ[key] = os.environ.get(key, "1")

    # Set cache path
    cache = Path(cache).expanduser()
    cache.mkdir(parents=True, exist_ok=True)

    # Get MSA server credentials from environment variables if not provided
    if use_msa_server:
        if msa_server_username is None:
            msa_server_username = os.environ.get("BOLTZ_MSA_USERNAME")
        if msa_server_password is None:
            msa_server_password = os.environ.get("BOLTZ_MSA_PASSWORD")
        if api_key_value is None:
            api_key_value = os.environ.get("MSA_API_KEY_VALUE")
        
        click.echo(f"MSA server enabled: {msa_server_url}")
        if api_key_value:
            click.echo("MSA server authentication: using API key header")
        elif msa_server_username and msa_server_password:
            click.echo("MSA server authentication: using basic auth")
        else:
            click.echo("MSA server authentication: no credentials provided")

    # Create output directories
    data = Path(data).expanduser()
    out_dir = Path(out_dir).expanduser()
    out_dir = out_dir / f"boltz_results_{data.stem}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Download necessary data and model
    if model == "boltz1":
        download_boltz1(cache)
    elif model == "boltz2":
        download_boltz2(cache)
    else:
        msg = f"Model {model} not supported. Supported: boltz1, boltz2."
        raise ValueError(f"Model {model} not supported.")

    # Validate inputs
    data = check_inputs(data)

    # Check method
    if method is not None:
        if model == "boltz1":
            msg = "Method conditioning is not supported for Boltz-1."
            raise ValueError(msg)
        if method.lower() not in const.method_types_ids:
            method_names = list(const.method_types_ids.keys())
            msg = f"Method {method} not supported. Supported: {method_names}"
            raise ValueError(msg)

    # Process inputs
    ccd_path = cache / "ccd.pkl"
    mol_dir = cache / "mols"
    process_inputs(
        data=data,
        out_dir=out_dir,
        ccd_path=ccd_path,
        mol_dir=mol_dir,
        use_msa_server=use_msa_server,
        msa_server_url=msa_server_url,
        msa_pairing_strategy=msa_pairing_strategy,
        msa_server_username=msa_server_username,
        msa_server_password=msa_server_password,
        api_key_header=api_key_header,
        api_key_value=api_key_value,
        boltz2=model == "boltz2",
        preprocessing_threads=preprocessing_threads,
        max_msa_seqs=max_msa_seqs,
    )

    # Load manifest
    manifest = Manifest.load(out_dir / "processed" / "manifest.json")

    # Filter out existing predictions
    filtered_manifest = filter_inputs_structure(
        manifest=manifest,
        outdir=out_dir,
        override=override,
    )

    # Load processed data
    processed_dir = out_dir / "processed"
    processed = BoltzProcessedInput(
        manifest=filtered_manifest,
        targets_dir=processed_dir / "structures",
        msa_dir=processed_dir / "msa",
        constraints_dir=(
            (processed_dir / "constraints")
            if (processed_dir / "constraints").exists()
            else None
        ),
        template_dir=(
            (processed_dir / "templates")
            if (processed_dir / "templates").exists()
            else None
        ),
        extra_mols_dir=(
            (processed_dir / "mols") if (processed_dir / "mols").exists() else None
        ),
    )

    # Set up trainer
    strategy = "auto"
    if (isinstance(devices, int) and devices > 1) or (
        isinstance(devices, list) and len(devices) > 1
    ):
        start_method = "fork" if platform.system() != "win32" and platform.system() != "Windows" else "spawn"
        strategy = DDPStrategy(start_method=start_method)
        if len(filtered_manifest.records) < devices:
            msg = (
                "Number of requested devices is greater "
                "than the number of predictions, taking the minimum."
            )
            click.echo(msg)
            if isinstance(devices, list):
                devices = devices[: max(1, len(filtered_manifest.records))]
            else:
                devices = max(1, min(len(filtered_manifest.records), devices))

    # Set up model parameters
    if model == "boltz2":
        diffusion_params = Boltz2DiffusionParams()
        step_scale = 1.5 if step_scale is None else step_scale
        diffusion_params.step_scale = step_scale
        pairformer_args = PairformerArgsV2()
    else:
        diffusion_params = BoltzDiffusionParams()
        step_scale = 1.638 if step_scale is None else step_scale
        diffusion_params.step_scale = step_scale
        pairformer_args = PairformerArgs()

    msa_args = MSAModuleArgs(
        subsample_msa=subsample_msa,
        num_subsampled_msa=num_subsampled_msa,
        use_paired_feature=model == "boltz2",
    )

    # Create prediction writer
    pred_writer = BoltzWriter(
        data_dir=processed.targets_dir,
        output_dir=out_dir / "predictions",
        output_format=output_format,
        boltz2=model == "boltz2",
        write_embeddings=write_embeddings,
    )

    # Set up trainer
    trainer = Trainer(
        default_root_dir=out_dir,
        strategy=strategy,
        callbacks=[pred_writer],
        accelerator=accelerator,
        devices=devices,
        precision=32 if model == "boltz1" else "bf16-mixed",
    )

    if filtered_manifest.records:
        msg = f"Running structure prediction for {len(filtered_manifest.records)} input"
        msg += "s." if len(filtered_manifest.records) > 1 else "."
        click.echo(msg)

        # Create data module
        if model == "boltz2":
            data_module = Boltz2InferenceDataModule(
                manifest=processed.manifest,
                target_dir=processed.targets_dir,
                msa_dir=processed.msa_dir,
                mol_dir=mol_dir,
                num_workers=num_workers,
                constraints_dir=processed.constraints_dir,
                template_dir=processed.template_dir,
                extra_mols_dir=processed.extra_mols_dir,
                override_method=method,
            )
        else:
            data_module = BoltzInferenceDataModule(
                manifest=processed.manifest,
                target_dir=processed.targets_dir,
                msa_dir=processed.msa_dir,
                num_workers=num_workers,
                constraints_dir=processed.constraints_dir,
            )

        # Load model
        if checkpoint is None:
            if model == "boltz2":
                checkpoint = cache / "boltz2_conf.ckpt"
            else:
                checkpoint = cache / "boltz1_conf.ckpt"

        predict_args = {
            "recycling_steps": recycling_steps,
            "sampling_steps": sampling_steps,
            "diffusion_samples": diffusion_samples,
            "max_parallel_samples": max_parallel_samples,
            "write_confidence_summary": True,
            "write_full_pae": write_full_pae,
            "write_full_pde": write_full_pde,
        }

        steering_args = BoltzSteeringParams()
        steering_args.fk_steering = use_potentials
        steering_args.physical_guidance_update = use_potentials
        steering_args.cdr3_steering = cdr3_steering
        steering_args.antigen_steering = antigen_steering
        steering_args.num_particles = num_particles

        # Validate CDR3 steering requires potentials
        if cdr3_steering and not use_potentials:
            msg = "CDR3 steering (--cdr3_steering) requires potentials to be enabled (--use_potentials)"
            raise click.UsageError(msg)

        # Validate antigen steering requires potentials
        if antigen_steering and not use_potentials:
            msg = "Antigen steering (--antigen_steering) requires potentials to be enabled (--use_potentials)"
            raise click.UsageError(msg)

        model_cls = Boltz2 if model == "boltz2" else Boltz1
        model_module = model_cls.load_from_checkpoint(
            checkpoint,
            strict=True,
            predict_args=predict_args,
            map_location="cpu",
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            use_kernels=not no_kernels,
            pairformer_args=asdict(pairformer_args),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
        )
        model_module.eval()

        # Run step-by-step mode if requested
        if step_by_step:
            if not _INSPECTION_AVAILABLE:
                msg = (
                    "Step-by-step mode requires the inspection module. "
                    "Make sure the inspection package is available."
                )
                raise click.ClickException(msg)

            click.echo("\nRunning step-by-step prediction mode\n")
            # Move model to GPU if available (model is loaded on CPU by default)
            if accelerator == "gpu" and torch.cuda.is_available():
                model_device = torch.device("cuda:0")
                model_module = model_module.to(model_device)
                click.echo(f"Model moved to {model_device}")
            else:
                model_device = torch.device("cpu")
            inspection_result = inspect_model_and_data(
                model_module, data_module, device=model_device
            )
            # Run submodules step by step for detailed inspection
            # Uses model's predict_args for recycling_steps and sampling_steps
            submodule_outputs = run_submodules_step_by_step(
                model=inspection_result['model'],
                batch=inspection_result['batch'],
                verbose=True,
            )
            # Write submodule outputs to from_submodules folder
            write_submodule_outputs(
                submodule_outputs=submodule_outputs,
                batch=inspection_result['batch'],
                output_dir=out_dir / "predictions",
                data_dir=processed.targets_dir,
                output_format=output_format,
                boltz2=(model == "boltz2"),
                write_embeddings=write_embeddings,
                write_confidence=True,
                write_pae=write_full_pae,
                write_pde=write_full_pde,
                verbose=True,
            )
        elif canonical_ensemble:
            import numpy as np
            import json
            from dataclasses import replace as dc_replace
            from boltz.data.types import Coords, Interface, StructureV2
            from boltz.data.write.mmcif import to_mmcif
            from boltz.data.write.pdb import to_pdb

            click.echo("\nRunning canonical ensemble sweep mode\n")

            # Move model to device
            if accelerator == "gpu" and torch.cuda.is_available():
                model_device = torch.device("cuda:0")
                model_module = model_module.to(model_device)
                click.echo(f"Model moved to {model_device}")
            else:
                model_device = torch.device("cpu")

            # Load batch from dataloader
            data_module.setup(stage="predict")
            dataloader = data_module.predict_dataloader()
            batch = next(iter(dataloader))
            batch = data_module.transfer_batch_to_device(batch, model_device, 0)

            # Extract canonical_ensemble_params from record
            record = batch["record"][0]
            ce_params = record.inference_options.canonical_ensemble_params
            if ce_params is None:
                msg = (
                    "--canonical_ensemble flag requires a canonical_ensemble "
                    "constraint in the YAML file."
                )
                raise click.ClickException(msg)

            beta_min, beta_max, beta_steps, top_k = ce_params
            beta_values = np.linspace(beta_min, beta_max, beta_steps)
            click.echo(
                f"Sweeping {beta_steps} beta values from {beta_min} to {beta_max}, "
                f"selecting top-{top_k}"
            )

            # ============================================================
            # TRUNK: Run once (input_embedder -> recycling -> MSA/pairformer)
            # ============================================================
            click.echo("Running trunk (input_embedder + recycling + MSA + pairformer)...")
            with torch.no_grad():
                s_inputs = model_module.input_embedder(batch)
                s_init = model_module.s_init(s_inputs)
                z_init = (
                    model_module.z_init_1(s_inputs)[:, :, None]
                    + model_module.z_init_2(s_inputs)[:, None, :]
                )
                relative_position_encoding = model_module.rel_pos(batch)
                z_init = z_init + relative_position_encoding
                z_init = z_init + model_module.token_bonds(batch["token_bonds"].float())
                if hasattr(model_module, 'bond_type_feature') and model_module.bond_type_feature:
                    z_init = z_init + model_module.token_bonds_type(batch["type_bonds"].long())
                z_init = z_init + model_module.contact_conditioning(batch)

                s = torch.zeros_like(s_init)
                z = torch.zeros_like(z_init)
                mask = batch["token_pad_mask"].float()
                pair_mask = mask[:, :, None] * mask[:, None, :]

                for i in range(recycling_steps + 1):
                    s = s_init + model_module.s_recycle(model_module.s_norm(s))
                    z = z_init + model_module.z_recycle(model_module.z_norm(z))
                    if hasattr(model_module, 'use_templates') and model_module.use_templates:
                        template_module = model_module.template_module
                        if hasattr(template_module, '_orig_mod'):
                            template_module = template_module._orig_mod
                        z = z + template_module(
                            z, batch, pair_mask, use_kernels=getattr(model_module, 'use_kernels', False)
                        )
                    msa_module = model_module.msa_module
                    if hasattr(msa_module, '_orig_mod'):
                        msa_module = msa_module._orig_mod
                    z = z + msa_module(
                        z, s_inputs, batch, use_kernels=getattr(model_module, 'use_kernels', False)
                    )
                    pairformer_module = model_module.pairformer_module
                    if hasattr(pairformer_module, '_orig_mod'):
                        pairformer_module = pairformer_module._orig_mod
                    s, z = pairformer_module(
                        s, z, mask=mask, pair_mask=pair_mask,
                        use_kernels=getattr(model_module, 'use_kernels', False),
                    )

                # Distogram (beta-independent)
                pdistogram = model_module.distogram_module(z)

            click.echo("Trunk complete.")

            # ============================================================
            # SWEEP: For each beta, run diffusion_conditioning + sampling + confidence
            # ============================================================
            all_results = []
            for beta_idx, beta_val in enumerate(beta_values):
                click.echo(f"  Beta {beta_idx+1}/{beta_steps}: beta={beta_val:.3f}")
                batch["cdr3_beta_value"] = torch.tensor(beta_val, device=model_device, dtype=torch.float32)

                with torch.no_grad():
                    # Diffusion conditioning (applies z *= (1+beta) on CDR regions)
                    q, c, to_keys, atom_enc_bias, atom_dec_bias, token_trans_bias = (
                        model_module.diffusion_conditioning(
                            s_trunk=s,
                            z_trunk=z,
                            relative_position_encoding=relative_position_encoding,
                            feats=batch,
                        )
                    )
                    diffusion_conditioning = {
                        'q': q, 'c': c, 'to_keys': to_keys,
                        'atom_enc_bias': atom_enc_bias,
                        'atom_dec_bias': atom_dec_bias,
                        'token_trans_bias': token_trans_bias,
                    }

                    # Structure module sampling
                    with torch.autocast("cuda", enabled=False):
                        struct_out = model_module.structure_module.sample(
                            s_trunk=s.float(),
                            s_inputs=s_inputs.float(),
                            feats=batch,
                            num_sampling_steps=sampling_steps,
                            atom_mask=batch["atom_pad_mask"].float(),
                            multiplicity=diffusion_samples,
                            max_parallel_samples=max_parallel_samples,
                            steering_args=getattr(model_module, 'steering_args', None),
                            diffusion_conditioning=diffusion_conditioning,
                        )

                    x_pred = struct_out['sample_atom_coords']

                    # Confidence module
                    confidence_out = model_module.confidence_module(
                        s_inputs=s_inputs.detach(),
                        s=s.detach(),
                        z=z.detach(),
                        x_pred=x_pred.detach(),
                        feats=batch,
                        pred_distogram_logits=pdistogram[:, :, :, 0].detach(),
                        multiplicity=diffusion_samples,
                        run_sequentially=True,
                        use_kernels=getattr(model_module, 'use_kernels', False),
                    )

                # Store results on CPU
                num_samples = x_pred.shape[0] if len(x_pred.shape) == 3 else 1
                for sample_idx in range(num_samples):
                    coords_sample = x_pred[sample_idx].cpu() if num_samples > 1 else x_pred.squeeze(0).cpu()

                    # Compute confidence score
                    complex_plddt = confidence_out.get('complex_plddt', torch.zeros(1))
                    iptm = confidence_out.get('iptm', torch.zeros(1))
                    ptm = confidence_out.get('ptm', torch.zeros(1))
                    if torch.is_tensor(iptm) and not torch.allclose(iptm, torch.zeros_like(iptm)):
                        tm_score = iptm
                    else:
                        tm_score = ptm

                    if complex_plddt.numel() > sample_idx:
                        c_plddt = complex_plddt[sample_idx].item() if complex_plddt.dim() > 0 else complex_plddt.item()
                    else:
                        c_plddt = complex_plddt.item()
                    if tm_score.numel() > sample_idx:
                        t_val = tm_score[sample_idx].item() if tm_score.dim() > 0 else tm_score.item()
                    else:
                        t_val = tm_score.item()

                    conf_score = (4 * c_plddt + t_val) / 5

                    # Extract per-sample confidence metrics
                    conf_dict = {'beta': float(beta_val), 'confidence_score': conf_score}
                    for key in ["ptm", "iptm", "ligand_iptm", "protein_iptm",
                                "complex_plddt", "complex_iplddt", "complex_pde", "complex_ipde"]:
                        if key in confidence_out:
                            val = confidence_out[key]
                            if torch.is_tensor(val):
                                if val.numel() > sample_idx:
                                    conf_dict[key] = val[sample_idx].item() if val.dim() > 0 else val.item()
                                else:
                                    conf_dict[key] = val.item()

                    plddt_sample = None
                    if 'plddt' in confidence_out:
                        plddt_tensor = confidence_out['plddt']
                        if plddt_tensor.dim() > 1 and plddt_tensor.shape[0] > sample_idx:
                            plddt_sample = plddt_tensor[sample_idx].cpu()
                        elif plddt_tensor.dim() == 1:
                            plddt_sample = plddt_tensor.cpu()

                    # Handle pair_chains_iptm
                    if "pair_chains_iptm" in confidence_out:
                        pair_iptm = confidence_out["pair_chains_iptm"]
                        conf_dict["chains_ptm"] = {
                            str(idx): pair_iptm[idx][idx][sample_idx].item()
                            if pair_iptm[idx][idx].numel() > sample_idx
                            else pair_iptm[idx][idx].item()
                            for idx in pair_iptm
                        }
                        conf_dict["pair_chains_iptm"] = {
                            str(idx1): {
                                str(idx2): pair_iptm[idx1][idx2][sample_idx].item()
                                if pair_iptm[idx1][idx2].numel() > sample_idx
                                else pair_iptm[idx1][idx2].item()
                                for idx2 in pair_iptm[idx1]
                            }
                            for idx1 in pair_iptm
                        }

                    all_results.append({
                        'coords': coords_sample,
                        'confidence_score': conf_score,
                        'confidence_dict': conf_dict,
                        'plddt': plddt_sample,
                        'beta': float(beta_val),
                    })

            # ============================================================
            # SELECT TOP-K and WRITE OUTPUT
            # ============================================================
            all_results.sort(key=lambda x: x['confidence_score'], reverse=True)
            selected = all_results[:top_k]

            click.echo(f"\nTop-{top_k} results:")
            for rank, res in enumerate(selected):
                click.echo(
                    f"  Rank {rank}: beta={res['beta']:.3f}, "
                    f"confidence={res['confidence_score']:.4f}"
                )

            # Write output files
            struct_dir = out_dir / "predictions" / record.id
            struct_dir.mkdir(parents=True, exist_ok=True)

            # Load structure for coordinate writing
            structure_path = processed.targets_dir / f"{record.id}.npz"
            structure: StructureV2 = StructureV2.load(structure_path)
            chain_map = {}
            for i, msk in enumerate(structure.mask):
                if msk:
                    chain_map[len(chain_map)] = i
            structure = structure.remove_invalid_chains()

            pad_mask = batch["atom_pad_mask"][0] if batch["atom_pad_mask"].dim() > 1 else batch["atom_pad_mask"]

            for rank, res in enumerate(selected):
                model_coord = res['coords']
                coord_unpad = model_coord[pad_mask.bool().cpu()]
                coord_unpad_np = coord_unpad.numpy()

                # Build structure
                atoms = structure.atoms.copy()
                atoms["coords"] = coord_unpad_np
                atoms["is_present"] = True
                coord_unpad_typed = [(x,) for x in coord_unpad_np]
                coord_unpad_typed = np.array(coord_unpad_typed, dtype=Coords)

                residues = structure.residues.copy()
                residues["is_present"] = True
                interfaces = np.array([], dtype=Interface)
                new_structure = dc_replace(
                    structure,
                    atoms=atoms,
                    residues=residues,
                    interfaces=interfaces,
                    coords=coord_unpad_typed,
                )

                outname = f"{record.id}_model_{rank}"
                plddts = res['plddt']

                if output_format == "pdb":
                    path = struct_dir / f"{outname}.pdb"
                    with path.open("w") as f:
                        f.write(to_pdb(new_structure, plddts=plddts, boltz2=(model == "boltz2")))
                else:
                    path = struct_dir / f"{outname}.cif"
                    with path.open("w") as f:
                        f.write(to_mmcif(new_structure, plddts=plddts, boltz2=(model == "boltz2")))

                # Write confidence JSON
                conf_path = struct_dir / f"confidence_{record.id}_model_{rank}.json"
                with conf_path.open("w") as f:
                    f.write(json.dumps(res['confidence_dict'], indent=4))

                # Write plddt
                if plddts is not None:
                    plddt_path = struct_dir / f"plddt_{record.id}_model_{rank}.npz"
                    plddt_np = plddts.numpy() if torch.is_tensor(plddts) else plddts
                    np.savez_compressed(plddt_path, plddt=plddt_np)

            click.echo(f"\nCanonical ensemble results written to {struct_dir}")
        else:
            # Compute structure predictions using normal forward pass
            trainer.predict(
                model_module,
                datamodule=data_module,
                return_predictions=False,
            )

    # Check if affinity predictions are needed
    if any(r.affinity for r in manifest.records):
        # Print header
        click.echo("\nPredicting property: affinity\n")

        # Validate inputs
        manifest_filtered = filter_inputs_affinity(
            manifest=manifest,
            outdir=out_dir,
            override=override,
        )
        if not manifest_filtered.records:
            click.echo("Found existing affinity predictions for all inputs, skipping.")
            return

        msg = f"Running affinity prediction for {len(manifest_filtered.records)} input"
        msg += "s." if len(manifest_filtered.records) > 1 else "."
        click.echo(msg)

        pred_writer = BoltzAffinityWriter(
            data_dir=processed.targets_dir,
            output_dir=out_dir / "predictions",
        )

        data_module = Boltz2InferenceDataModule(
            manifest=manifest_filtered,
            target_dir=out_dir / "predictions",
            msa_dir=processed.msa_dir,
            mol_dir=mol_dir,
            num_workers=num_workers,
            constraints_dir=processed.constraints_dir,
            template_dir=processed.template_dir,
            extra_mols_dir=processed.extra_mols_dir,
            override_method="other",
            affinity=True,
        )

        predict_affinity_args = {
            "recycling_steps": 5,
            "sampling_steps": sampling_steps_affinity,
            "diffusion_samples": diffusion_samples_affinity,
            "max_parallel_samples": 1,
            "write_confidence_summary": False,
            "write_full_pae": False,
            "write_full_pde": False,
        }

        # Load affinity model
        if affinity_checkpoint is None:
            affinity_checkpoint = cache / "boltz2_aff.ckpt"

        steering_args = BoltzSteeringParams()
        steering_args.fk_steering = False
        steering_args.physical_guidance_update = False
        steering_args.contact_guidance_update = False
        
        model_module = Boltz2.load_from_checkpoint(
            affinity_checkpoint,
            strict=True,
            predict_args=predict_affinity_args,
            map_location="cpu",
            diffusion_process_args=asdict(diffusion_params),
            ema=False,
            pairformer_args=asdict(pairformer_args),
            msa_args=asdict(msa_args),
            steering_args=asdict(steering_args),
            affinity_mw_correction=affinity_mw_correction,
        )
        model_module.eval()

        trainer.callbacks[0] = pred_writer
        trainer.predict(
            model_module,
            datamodule=data_module,
            return_predictions=False,
        )


if __name__ == "__main__":
    cli()
