# contact_constrained trainer
boltz predict examples/contact_constrained/7Y0O_HLA.yaml --output_format pdb

# unconstrained trainer
boltz predict examples/unconstrained/7Y0O_HLA.yaml --output_format pdb

# unconstrained step_by_step
boltz predict examples/unconstrained/7Y0O_HLA.yaml --output_format pdb --step_by_step

# contact_constrained step_by_step
boltz predict examples/contact_constrained/7Y0O_HLA.yaml --output_format pdb --step_by_step

# steering_cdr3 step_by_step
boltz predict examples/contact_constrained/7Y0O_HLA_steer_cdr3.yaml --output_format pdb --use_potentials --cdr3_steering --step_by_step --diffusion_samples 5

# steering_cdr3 trainer
boltz predict examples/contact_constrained/7Y0O_HLA_steer_cdr3.yaml --output_format pdb --use_potentials --cdr3_steering --diffusion_samples 5

# antigen_steering
boltz predict examples/contact_constrained/7Y0O_HLA_antigen_steer.yaml --output_format pdb --diffusion_samples 5 --use_potentials  --antigen_steering