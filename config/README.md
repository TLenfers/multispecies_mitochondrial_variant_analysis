# General Configuration
To configure this workflow, modify [`config.yaml`](/config/config.yaml) according to your needs, following the explanations provided in the file.

## Sample and unit sheet
Add samples to [`config.yaml`](/config/config.yaml) . For each sample, the columns sample_name has to be defined.

## Reference
Specify the reference genome to use. Is dog, human, or mouse specified, the corresponding reference genome will be loaded.
Otherwise define the name of the reference in the config file and put the reference file in the `workflow/data/reference`.
