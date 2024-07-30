#! /usr/bin/env python3
from orthologue_analysis import (
    parse_args,
    main
)
from orthologue_analysis.species import (
    SpeciesList,
    HaematobiumClade,
    MansoniClade,
    JaponicumClade,
    IndicumClade,
    NewSpeciesClade,
)


if __name__ == "__main__":
    args = parse_args()
    
    SPECIES_LIST = SpeciesList([
        HaematobiumClade("bovis", "TD2_PRJEB44434"),
        HaematobiumClade("curassoni", "PRJEB44434"),
        HaematobiumClade("guineensis", "PRJEB44434"),
        HaematobiumClade("haematobium", "TD2_PRJEB44434"),
        HaematobiumClade("intercalatum", "TD2_PRJEB44434"),
        HaematobiumClade("margrebowiei", "PRJEB44434"),
        HaematobiumClade("mattheei", "PRJEB44434"),
        MansoniClade("rodhaini", "TD2_PRJEB44434"),
        MansoniClade("mansoni", "PRJEA36577"),
        JaponicumClade("japonicum", "PRJNA520774"),
        IndicumClade("spindale", "PRJEB44434"),
        NewSpeciesClade("turkestanicum", "PRJEB44434")
        ],
        **vars(args)
    )

    main(args, SPECIES_LIST)
