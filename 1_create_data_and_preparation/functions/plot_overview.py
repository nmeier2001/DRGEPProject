import matplotlib.pyplot as plt

class PlotOverview:
    '''
    Creates an overview of the study results, including plots for conversion vs. SEI, energy efficiency vs. conversion, SEI vs. energy efficiency, electron density vs. time, conversion vs. time, TGAS vs. time, CO2 vs. time, CO vs. time, and ratios of different species to CO2 vs. time.

    :param study: dictionary containing the study data
    :param path_save: path to save the plots
    :param conversion_list: list of conversion values for each run
    :param SEI_list: list of SEI values for each run
    :param energy_efficiency_list: list of energy efficiency values for each run
    :param species: list of species to plot ratios against CO2
    :param is_rate_study: boolean indicating if the study is a Combined Rate Study or not
    '''
    def __init__(self, study, path_save: str = None, conversion_list: list = None, SEI_list: list = None, energy_efficiency_list: list = None, species: list = None, is_rate_study: bool = False):
        self.study = study
        self.path_save = path_save
        self.conversion_list = conversion_list
        self.SEI_list = SEI_list
        self.energy_efficiency_list = energy_efficiency_list
        self.species = species
        self.is_rate_study = is_rate_study
        plt.rcParams.update({'font.size': 14})
        self.plot_conversion_vs_SEI()
        self.plot_energy_efficiency_vs_conversion()
        self.plot_SEI_vs_energy_efficiency()
        self.plot_e_density_vs_time()
        self.plot_conversion_vs_time()
        self.plot_TGAS_vs_time()
        self.plot_CO2_vs_time()
        self.plot_CO_vs_time()
        self.plot_species_to_CO2_ratio_vs_time(species)
        self.plot_EN_vs_time()
        self.plot_total_CO2_density_vs_time()
        

    # plot conversion vs. SEI for each run
    def plot_conversion_vs_SEI(self):
        SEI_list = self.SEI_list
        conversion_list = self.conversion_list
        plt.figure(figsize=(8, 6))
        plt.scatter(SEI_list, conversion_list, alpha=0.7)
        plt.xlabel('SEI [eV]')
        plt.ylabel('Konversion CO2')
        plt.title('Konversion vs. SEI')
        plt.xscale('log')
        #plt.yscale('log')
        plt.grid(True)
        plt.savefig(self.path_save / f'Conversion_vs_SEI.png')

    #plot energy efficiency vs. conversion for each run
    def plot_energy_efficiency_vs_conversion(self):
        plt.figure(figsize=(8, 6))
        plt.scatter(self.energy_efficiency_list, self.conversion_list, alpha=0.7)
        plt.xlabel('Energieeffizienz')
        plt.ylabel('Konversion')
        plt.title('Konversion vs. Energieeffizienz')
        #plt.xscale('log')
        #plt.yscale('log')
        plt.grid(True)
        plt.savefig(self.path_save / f'Conversion_vs_Energieeffizienz.png')

    # plot SEI vs. energy efficiency
    def plot_SEI_vs_energy_efficiency(self):
        plt.figure(figsize=(8, 6))
        plt.scatter(self.SEI_list, self.energy_efficiency_list, alpha=0.7)
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel('SEI [eV]')
        plt.ylabel('Energieeffizienz')
        plt.title('SEI vs. Energieeffizienz')
        plt.grid(True)
        plt.savefig(self.path_save / f'SEI_vs_Energieeffizienz.png')

    # plot e_density vs. time for each run in the same figure
    def plot_e_density_vs_time(self):
            plt.figure(figsize=(8, 6))
            for idx, run in self.study['runs'].items():
                df = run['tec']
                plt.plot(df['T (s)'], df['E'], label=f"{run['param_value']}")
            plt.xscale('log')
            plt.yscale('log')
            plt.xlabel('Zeit [s]')
            plt.ylabel('Elektronendichte')
            # plt.title('Elektronendichte vs. Zeit')
            plt.legend(title = "Plug Power [W/cm³]")
            plt.grid(True)
            plt.savefig(self.path_save / f'Electron_Density_vs_Time.png')

    # plot conversion vs. time for each run in the same figure
    def plot_conversion_vs_time(self):
        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['CONVERSION'], label=f"{run['param_value']}")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('Konversion CO2')
        # plt.title('Konversion vs. Zeit')
        plt.legend(title = "Plug Power [W/cm³]")
        plt.grid(True)
        plt.savefig(self.path_save / f'Conversion_vs_Time.png')

    # plot TGAS vs. time for each run in the same figure
    def plot_TGAS_vs_time(self):
        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['TGAS'], label=f"{run['param_value']}")
        plt.xscale('log')
        #plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('TGAS')
        # plt.title('TGAS vs. Zeit')
        plt.legend(title = "Leistung [W]")
        plt.grid(True)
        plt.savefig(self.path_save / f'TGAS_vs_Time.png')

    # plot CO vs. time for each run in the same figure
    def plot_CO_vs_time(self):
        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['CO'], label=f"{run['param_value']}")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('CO [m⁻³]')
        # plt.title('CO vs. Zeit')
        plt.legend(title = "TCYCLE")
        plt.grid(True)
        plt.savefig(self.path_save / f'CO_vs_Time.png')

    # plot ratios of different species to CO2 vs. time for each run in the same figure
    def plot_species_to_CO2_ratio_vs_time(self, species: list):
        plt.figure(figsize=(8, 6))
        for molecule in species:
            plt.plot(self.study['runs'][1]['tec']['T (s)'], self.study['runs'][1]['tec'][f'{molecule}_to_CO2_ratio'], label=f"{molecule}")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('Verhältnis Molekül/CO2')
        # plt.title('Verhältnis Molekül/CO2 vs. Zeit')
        plt.legend(title = "Molekül")
        plt.grid(True)
        plt.savefig(self.path_save / f'Molecule_to_CO2_ratio_vs_Time_run6.png')

        plt.figure(figsize=(8, 6))

    # plot CO2 vs time for each run in the same figure
    def plot_CO2_vs_time(self):
        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['CO2'], label=f"{run['param_value']}")
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('CO2 [m⁻³]')
        # plt.title('CO2 vs. Zeit')
        plt.legend(title = "Leistung [W]")
        plt.grid(True)
        plt.savefig(self.path_save / 'CO2_vs_Time.png', dpi=300,  bbox_inches='tight')

    def plot_total_CO2_density_vs_time(self, co2_species: list = None):
        plt.figure(figsize=(11, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
 
            if co2_species is not None:
                # explizit übergebene Spaltennamen verwenden
                cols = [c for c in co2_species if c in df.columns]
            else:
                # automatische Erkennung: alle Spalten, die mit 'CO2' beginnen
                # und keine Verhältnis-Spalten (_to_CO2_ratio) sind
                cols = [c for c in df.columns
                        if c.startswith('CO2') and not c.endswith('_ratio')]
 
            total_CO2 = df[cols].sum(axis=1)
            plt.plot(df['T (s)'], total_CO2, label=f"{run['param_value']}")
 
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('CO2-Dichte [m⁻³]')
        # plt.title('Gesamte CO2-Dichte vs. Zeit')
        plt.legend(title="Leistung [W]")
        plt.grid(True)
        plt.savefig(self.path_save / 'Total_CO2_Density_vs_Time.png', dpi=400,  bbox_inches='tight')

    def plot_EN_vs_time(self):
        # ── EN-Kurven ─────────────────────────────────────────────────────────────────
        fig, ax = plt.subplots(figsize=(10, 5))

        for idx, run in self.study['runs'].items():
            #print(f"Run {idx}: param_value = {run['param_value']}")
            df = run['tec']
            ax.plot(df['T (s)'], df['EN_Td'], label=f"{run['param_value']}")

        # ── Puls aus erstem Run (Nam-Parameter) ───────────────────────────────────────
        nam     = self.study['runs'][1]['nam']
        powplug = nam.get('POWPLUG')  # [0.0, 2e7, 0.0, 0.0]
        xpow    = nam.get('XPOW')     # [0., 15e-9, 30e-9, 2e-4, tcycle]

        if powplug and xpow and isinstance(powplug, list) and isinstance(xpow, list):
            ax2 = ax.twinx()

            # Dreieckspuls: x0 → peak → end
            x0    = xpow[0]
            peak  = xpow[1]
            end   = xpow[2]
            power = powplug[1]  # Peak-Leistung

            t = [x0,    peak,  end,  xpow[3]]
            p = [0.0,   power, 0.0,  0.0    ]

            ax2.plot(t, p, color='gray', linestyle='--', linewidth=1.2, label='Puls')
            ax2.set_ylabel('EN [Td]')
            ax2.set_ylim(0, power * 3)
            # ax2.set_xlim(0, end)
            ax2.legend(loc='lower right')

        ax.set_xscale('log')
        ax.set_xlabel('Zeit [s]')
        ax.set_ylabel('EN [Td]')
        ax.set_title('EN vs. Zeit')
        ax.legend(title="EN [Td]", loc='upper right')
        ax.grid(True)
        plt.tight_layout()
        plt.savefig(self.path_save / 'EN_vs_Time.png', dpi=300)


        # Plot TEC data for EN_Td vs. time for each run in the same figure
        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['TE'], label=f"{run['param_value']}")
        plt.xscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('TE [eV]')
        # plt.title('TE vs. Zeit')
        plt.legend(title="TE [K]")
        plt.grid(True)
        plt.savefig(self.path_save / 'TE_vs_Time.png', dpi=400)

        plt.figure(figsize=(8, 6))
        for idx, run in self.study['runs'].items():
            df = run['tec']
            plt.plot(df['T (s)'], df['E'], label=f"{run['param_value']}")
        plt.xscale('log')
        plt.xlabel('Zeit [s]')
        plt.ylabel('E [cm⁻³]')
        # plt.title('E vs. Zeit')
        plt.legend(title="E [cm⁻³]")
        plt.grid(True)
        plt.savefig(self.path_save / 'E_vs_Time.png', dpi=400)