#! /usr/bin/env python3
import os.path

from orthologue_analysis import (
    parse_args,
    main
)
from orthologue_analysis.species import (
    SpeciesList,
    HaematobiumCladeFromTool,
    MansoniCladeFromTool,
    JaponicumCladeFromTool,
)

if __name__ == "__main__":
    mars_data_dir = os.path.join("data", "from_MARS", "")
    args = parse_args()
    args.prefix_cut = "transcript_"

    SPECIES_LIST = SpeciesList([
        HaematobiumCladeFromTool("bovis", data_dir=mars_data_dir, data_label="Sbov_LT", prot_filename_suffix=".fa"),
        HaematobiumCladeFromTool("curassoni", data_dir=mars_data_dir, data_label="Scur_LT", prot_filename_suffix=".fa"),
        HaematobiumCladeFromTool("guineensis", data_dir=mars_data_dir, data_label="Sgui_LT", prot_filename_suffix=".fa"),
        HaematobiumCladeFromTool("haematobium", data_dir=mars_data_dir, data_label="Shae_LT", prot_filename_suffix=".fa"),
        HaematobiumCladeFromTool("intercalatum", data_dir=mars_data_dir, data_label="Sint_LT", prot_filename_suffix=".fa"),
        MansoniCladeFromTool("rodhaini", data_dir=mars_data_dir, data_label="Srod_LT", prot_filename_suffix=".fa"),
        MansoniCladeFromTool("mansoni", data_dir=mars_data_dir, data_label="Sman_LT", prot_filename_suffix=".fa"),
        MansoniCladeFromTool("mansoni_braker3_reann", data_dir=mars_data_dir, data_label="Sman_braker3_LT", prot_filename_suffix=".fa", skip_plot=True),
        JaponicumCladeFromTool("japonicum", data_dir=mars_data_dir, data_label="Sjap_LT", prot_filename_suffix=".fa"),
        ],
        **vars(args)
    )

    main(args, SPECIES_LIST)
