'''
info_reduced = {
    powplug: {
        epsilon: {
            kept_species: [list of species],
            kept_reactions: number of reactions,
            num_species: number of species,
    },
}
}
'''

def create_data_reduced_dic(kept_species, eliminated_species, number_reactions, powplug, epsilon):
    info_reduced = {}
    info_reduced['kept_species'] = kept_species
    info_reduced['eliminated_species'] = eliminated_species
    info_reduced['kept_reactions'] = number_reactions
    info_reduced['num_species'] = len(kept_species)

    return info_reduced