import json
from pathlib import Path


def load_knowledge_base(path=None):
    dataset_path = Path(path) if path else Path(__file__).with_name("datasets.json")
    with dataset_path.open(encoding="utf-8") as dataset_file:
        return json.load(dataset_file)


def describe_smiles(smiles):
    from rdkit import Chem
    from rdkit.Chem import rdMolDescriptors

    molecule = Chem.MolFromSmiles(smiles.strip())
    if molecule is None:
        raise ValueError("La formule SMILES est invalide.")

    atoms = list(molecule.GetAtoms())
    ring_nitrogens = sum(
        atom.GetSymbol() == "N" and atom.IsInRing()
        for atom in atoms
    )
    halogens = {
        symbol: sum(atom.GetSymbol() == symbol for atom in atoms)
        for symbol in ("F", "Cl", "Br", "I")
    }
    detected_features = [
        name for name, pattern in FUNCTIONAL_GROUP_SMARTS.items()
        if molecule.HasSubstructMatch(Chem.MolFromSmarts(pattern))
    ]

    return {
        "formula": rdMolDescriptors.CalcMolFormula(molecule),
        "molecular_weight": rdMolDescriptors.CalcExactMolWt(molecule),
        "clogp": rdMolDescriptors.CalcCrippenDescriptors(molecule)[0],
        "tpsa": rdMolDescriptors.CalcTPSA(molecule),
        "hbd": rdMolDescriptors.CalcNumHBD(molecule),
        "hba": rdMolDescriptors.CalcNumHBA(molecule),
        "rotatable_bonds": rdMolDescriptors.CalcNumRotatableBonds(molecule),
        "ring_count": rdMolDescriptors.CalcNumRings(molecule),
        "aromatic_ring_count": rdMolDescriptors.CalcNumAromaticRings(molecule),
        "aromatic": any(atom.GetIsAromatic() for atom in atoms),
        "nitrogen_count": sum(atom.GetSymbol() == "N" for atom in atoms),
        "ring_nitrogen_count": ring_nitrogens,
        "oxygen_count": sum(atom.GetSymbol() == "O" for atom in atoms),
        "halogens": halogens,
        "detected_features": detected_features,
    }


FEATURE_LABELS = {
    "core_scaffold": "Squelette central",
    "nitrogen_count": "Azotes cycliques",
    "aromaticity": "Aromaticite",
    "aromatic_surface": "Cycles aromatiques",
    "lipophilicity": "Lipophilie",
    "polarity": "Polarite",
    "hydrogen_bond_donors": "Donneurs de liaison hydrogene",
    "hydrogen_bond_acceptors": "Accepteurs de liaison hydrogene",
    "molecular_flexibility": "Flexibilite",
    "amide_linker": "Liaison amide detectee",
    "ester_linker": "Liaison ester detectee",
    "carbonyl_group": "Groupe carbonyle detecte",
    "halogenation": "Halogenes",
    "fluorination": "Fluoration",
}

EFFECT_LABELS = {
    "changes_aromatic_geometry": "geometrie aromatique",
    "changes_electronic_distribution": "distribution electronique",
    "changes_hydrogen_bonding_pattern": "liaisons hydrogene",
    "changes_receptor_interactions": "interactions avec le recepteur",
    "electronic_properties": "proprietes electroniques",
    "polarity": "polarite",
    "hydrogen_bonding": "liaisons hydrogene",
    "binding_geometry": "geometrie de liaison",
    "pi_interactions": "interactions pi",
    "hydrophobic_interactions": "interactions hydrophobes",
    "rigidity": "rigidite",
    "conformational_constraint": "contrainte conformationnelle",
    "hydrophobic_contact": "contacts hydrophobes",
    "shape_complementarity": "complementarite de forme",
    "membrane_partitioning": "repartition membranaire",
    "distribution": "distribution",
    "protein_binding": "liaison aux proteines",
    "solubility": "solubilite",
    "membrane_permeability": "permeabilite membranaire",
    "receptor_interactions": "interactions avec le recepteur",
    "receptor_hydrogen_bonding": "liaisons hydrogene au recepteur",
    "conformational_entropy": "entropie conformationnelle",
    "binding_orientation": "orientation de liaison",
    "receptor_fit": "ajustement au recepteur",
    "electronic_effects": "effets electroniques",
    "lipophilicity": "lipophilie",
    "metabolic_stability": "stabilite metabolique",
    "polarizability": "polarisabilite",
    "molecular_size": "taille moleculaire",
    "ionization": "ionisation",
    "hydrogen_bond_donors": "donneurs de liaison hydrogene",
    "hydrogen_bond_acceptors": "accepteurs de liaison hydrogene",
    "ester_groups": "groupes ester",
    "amide_groups": "groupes amide",
    "benzylic_positions": "positions benzyliques",
    "heteroatoms": "heteroatomes",
    "free_fraction": "fraction libre",
    "duration": "duree d'action",
}

SCOPE_LABELS = {
    "family_dependent": "depend du squelette chimique",
    "scaffold_dependent": "depend du squelette chimique",
    "regioisomer_dependent": "depend de la position des substituants",
    "structure_dependent": "depend de la structure",
    "position_and_geometry_dependent": "depend de la position et de la geometrie",
    "series_dependent": "depend de la serie et de l'essai",
    "specific_series": "limitee a certaines series",
    "specific_scaffolds_and_assays": "limitee a certains squelettes et essais",
    "substituent_dependent": "depend du substituant",
}

CAVEAT_LABELS = {
    "larger_aromatic_surface_is_not_automatically_better":
        "une surface aromatique plus grande n'est pas automatiquement meilleure",
    "binding_affinity_and_membrane_partitioning_are_distinct_properties":
        "affinite et repartition membranaire sont des proprietes distinctes",
    "lipophilicity_alone_does_not_predict_CB1_affinity":
        "la lipophilie seule ne predit pas l'affinite CB1",
    "assay_dependent": "depend du systeme experimental",
    "assay_and_system_dependent": "depend de l'essai et du systeme biologique",
    "ionization_state_depends_on_pH": "l'etat d'ionisation depend du pH",
    "effect_depends_on_position_and_scaffold":
        "l'effet depend de la position et du squelette",
    "fluorination_does_not_automatically_increase_affinity_or_potency":
        "la fluoration n'augmente pas automatiquement l'affinite ou la puissance",
    "the_relationship_is_not_a_universal_4_to_6_rule":
        "il n'existe pas de regle universelle de longueur 4 a 6 carbones",
    "stereoisomers_can_have_different_pharmacology":
        "les stereoisomeres peuvent avoir des profils pharmacologiques differents",
}

RULE_LABELS = {
    "rule_01": "Ne pas deduire l'affinite CB1 de la seule lipophilie.",
    "rule_02": "Ne pas deduire EC50 a partir de Ki seul.",
    "rule_03": "Ne pas deduire Emax a partir de Ki seul.",
    "rule_04": "Ne pas generaliser une observation entre series chimiques differentes.",
    "rule_05": "Distinguer pharmacologie du recepteur et pharmacocinetique.",
    "rule_06": "Distinguer l'activite de la molecule mere de celle des metabolites.",
    "rule_07": "Tenir compte du systeme experimental pour comparer EC50 ou Emax.",
    "rule_08": "Considerer stereochimie et regiochemie.",
    "rule_09": "Une molecule plus lipophile n'est pas automatiquement un agoniste CB1 plus puissant.",
    "rule_10": "Un changement structural peut affecter plusieurs proprietes a la fois.",
}

FUNCTIONAL_GROUP_SMARTS = {
    "amide_linker": "[CX3](=[OX1])[NX3]",
    "ester_linker": "[CX3](=[OX1])[OX2][#6]",
    "carbonyl_group": "[CX3]=[OX1]",
}

METRIC_LABELS = {
    "binding_affinity": "affinite de liaison",
    "functional_potency": "puissance fonctionnelle",
    "maximum_functional_response": "reponse fonctionnelle maximale",
}

METRIC_INTERPRETATIONS = {
    "lower_value_generally_means_higher_affinity": "une valeur plus faible indique en general une affinite plus elevee",
    "lower_value_generally_means_higher_potency": "une valeur plus faible indique en general une puissance plus elevee",
    "describes_efficacy": "decrit l'efficacite maximale, pas l'affinite",
}

RELATIONSHIP_LABELS = {
    "affinity_vs_potency": "Ki et EC50 mesurent des proprietes differentes.",
    "affinity_vs_efficacy": "Une forte affinite ne signifie pas necessairement une forte efficacite.",
    "potency_vs_efficacy": "EC50 et Emax sont deux dimensions independantes.",
    "full_agonism": "La reponse maximale d'un agoniste complet depend de l'essai compare.",
    "biased_signaling": "Des ligands peuvent differencier les voies de signalisation activees.",
    "receptor_selectivity": "L'affinite CB1 seule ne determine pas la selectivite entre recepteurs.",
}

OBSERVATION_LABELS = {
    "indazole_vs_indole": "Certaines series derivees d'acides amines rapportent une activite CB1 differente entre indazole et indole.",
    "indole_chain_length": "Des etudes de series indole rapportent une dependance a la longueur de chaine hydrophobe.",
    "substituent_effects": "De petits changements structuraux peuvent modifier plusieurs mesures experimentales.",
    "fluorination_effects": "La fluoration peut modifier electronique, lipophilie et metabolisme.",
}

PROPERTY_LABELS = {
    "BBB_permeability": "Passage de la barriere hematoencephalique",
    "membrane_partitioning": "Repartition membranaire",
    "protein_binding": "Liaison aux proteines",
    "metabolic_stability": "Stabilite metabolique",
    "active_metabolites": "Metabolites actifs",
}


def build_knowledge_report(descriptors, knowledge):
    features = knowledge["structural_features"]
    lines = [
        "Profil calcule depuis le SMILES",
        f"Formule : {descriptors['formula']}",
        f"Masse molaire : {descriptors['molecular_weight']:.2f} g/mol",
        f"cLogP estime : {descriptors['clogp']:.2f}",
        f"TPSA : {descriptors['tpsa']:.2f} A^2",
        f"Donneurs / accepteurs de liaison H : {descriptors['hbd']} / {descriptors['hba']}",
        f"Cycles / cycles aromatiques : {descriptors['ring_count']} / {descriptors['aromatic_ring_count']}",
        f"Azotes totaux / cycliques : {descriptors['nitrogen_count']} / {descriptors['ring_nitrogen_count']}",
        f"Oxygenes : {descriptors['oxygen_count']}",
        f"Liaisons rotatives : {descriptors['rotatable_bonds']}",
        "Halogenes : " + (
            ", ".join(f"{element} x{count}" for element, count in descriptors["halogens"].items() if count)
            or "aucun"
        ),
        "",
        "Mise en contexte depuis datasets.json",
    ]

    feature_values = {
        "core_scaffold": ", ".join(features["core_scaffold"]["examples"]),
        "nitrogen_count": f"{descriptors['ring_nitrogen_count']} azote(s) cyclique(s)",
        "aromaticity": "presente" if descriptors["aromatic"] else "absente",
        "aromatic_surface": f"{descriptors['aromatic_ring_count']} cycle(s) aromatique(s)",
        "lipophilicity": f"cLogP estime = {descriptors['clogp']:.2f}",
        "polarity": f"TPSA {descriptors['tpsa']:.2f} A^2 ; HBD {descriptors['hbd']} ; HBA {descriptors['hba']}",
        "hydrogen_bond_donors": str(descriptors["hbd"]),
        "hydrogen_bond_acceptors": str(descriptors["hba"]),
        "molecular_flexibility": f"{descriptors['rotatable_bonds']} liaison(s) rotative(s)",
        "amide_linker": "present" if "amide_linker" in descriptors["detected_features"] else "absent",
        "ester_linker": "present" if "ester_linker" in descriptors["detected_features"] else "absent",
        "carbonyl_group": "present" if "carbonyl_group" in descriptors["detected_features"] else "absent",
        "halogenation": ", ".join(
            f"{element} x{count}" for element, count in descriptors["halogens"].items() if count
        ) or "aucun halogene detecte",
        "fluorination": f"{descriptors['halogens']['F']} atome(s) de fluor",
    }
    relevant_features = [
        "core_scaffold", "nitrogen_count", "aromaticity", "aromatic_surface",
        "lipophilicity", "polarity", "hydrogen_bond_donors",
        "hydrogen_bond_acceptors", "molecular_flexibility", "amide_linker",
        "ester_linker", "carbonyl_group",
    ]
    if any(descriptors["halogens"].values()):
        relevant_features.append("halogenation")
        if descriptors["halogens"]["F"]:
            relevant_features.append("fluorination")

    for key in relevant_features:
        feature = features.get(key)
        if not feature:
            continue
        lines.append(f"- {FEATURE_LABELS[key]} : {feature_values[key]}")
        if key in FUNCTIONAL_GROUP_SMARTS and key not in descriptors["detected_features"]:
            continue
        effects = feature.get("possible_effects", feature.get("effects", []))
        if effects:
            translated = [
                EFFECT_LABELS.get(effect) or effect.replace("_", " ")
                for effect in effects
            ]
            lines.append("  Facteurs evoques : " + ", ".join(translated))
        scope = feature.get("scope")
        if scope:
            lines.append("  Portee : " + (SCOPE_LABELS.get(scope) or scope.replace("_", " ")))
        caveat = feature.get("caveat")
        if caveat:
            lines.append("  Limite : " + (CAVEAT_LABELS.get(caveat) or caveat.replace("_", " ")))

    target = knowledge["targets"]["CB1"]
    parameter_labels = {
        "Ki": "affinite de liaison Ki",
        "EC50": "puissance fonctionnelle EC50",
        "Emax": "efficacite maximale Emax",
        "binding_mode": "mode de liaison",
        "functional_activity": "activite fonctionnelle",
    }
    lines.extend([
        "",
        f"Cible de reference : CB1 ({target['type']}, signalisation {'/'.join(target['signaling'])}).",
        "Parametres experimentaux distincts : " + ", ".join(
            parameter_labels.get(parameter) or parameter
            for parameter in target["main_parameters"]
        ) + ".",
        "Aucune valeur Ki, EC50 ou Emax n'est predite : le JSON ne contient pas de mesures associees a des molecules pour entrainer un modele.",
        "",
        "Regles d'interpretation du dataset",
    ])

    for key in knowledge["interpretation_rules"]:
        lines.append("- " + (RULE_LABELS.get(key) or key.replace("_", " ")))

    lines.extend(["", "Sens des mesures experimentales (dataset)"])
    for name, metric in knowledge["metrics"].items():
        if "meaning" not in metric:
            if name == "efficacy_class":
                lines.append("- Classes d'activite : " + ", ".join(metric["values"]))
            continue
        meaning = METRIC_LABELS.get(metric["meaning"]) or metric["meaning"].replace("_", " ")
        interpretation = METRIC_INTERPRETATIONS.get(metric.get("interpretation", ""))
        line = f"- {name} : {meaning}"
        if interpretation:
            line += f" ; {interpretation}"
        caveat = metric.get("caveat")
        if caveat:
            line += " ; " + (CAVEAT_LABELS.get(caveat) or caveat.replace("_", " "))
        lines.append(line)

    lines.extend(["", "Relations pharmacologiques du dataset"])
    for key, relationship in knowledge["pharmacology_relationships"].items():
        statement = RELATIONSHIP_LABELS.get(key) or key.replace("_", " ")
        caveat = relationship.get("caveat")
        if caveat:
            statement += " (" + (CAVEAT_LABELS.get(caveat) or caveat.replace("_", " ")) + ")"
        lines.append("- " + statement)

    lines.extend(["", "Observations de series (non universelles)"])
    for key, observation in knowledge["known_series_level_observations"].items():
        statement = OBSERVATION_LABELS.get(key) or key.replace("_", " ")
        scope = observation.get("scope")
        if scope:
            statement += " Portee : " + (SCOPE_LABELS.get(scope) or scope.replace("_", " ")) + "."
        lines.append("- " + statement)

    lines.extend([
        "",
        "Facteurs pharmacocinetiques (non predits depuis le SMILES)",
    ])
    for key, property_data in knowledge["pharmacokinetic_properties"].items():
        label = PROPERTY_LABELS.get(key) or key.replace("_", " ")
        influencers = property_data.get("influenced_by", [])
        line = "- " + label
        if influencers:
            translated = [EFFECT_LABELS.get(item) or item.replace("_", " ") for item in influencers]
            line += " ; facteurs cites : " + ", ".join(translated)
        caveat = property_data.get("caveat")
        if caveat:
            line += " ; " + (CAVEAT_LABELS.get(caveat) or caveat.replace("_", " "))
        lines.append(line)

    lines.append("")
    lines.append("Avertissement : observations educatives, dependantes des series et des essais ; pas des regles d'optimisation universelles.")
    return "\n".join(lines)