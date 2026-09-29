import math
import tkinter as tk
from tkinter import messagebox, ttk

from sar_engine import describe_smiles


class SimpleMoleculeTester(tk.Toplevel):
    STRUCTURE_FIELDS = [
        ("name", "Nom / identifiant", "text"),
        ("smiles", "SMILES (facultatif)", "text"),
        ("core_scaffold", "Squelette central", "scaffold"),
        ("substituents", "Substituants", "text"),
        ("linker", "Liaison / linker", "text"),
        ("nitrogen_count", "Azotes totaux", "count"),
        ("ring_nitrogen_count", "Azotes cycliques", "count"),
        ("heteroatom_position", "Position des heteroatomes", "text"),
        ("aromaticity", "Aromaticite", "yes_no"),
        ("ring_count", "Nombre de cycles", "count"),
        ("aromatic_surface", "Cycles aromatiques", "count"),
        ("hydrophobic_domain", "Domaine hydrophobe", "text"),
        ("molecular_weight", "Masse molaire (g/mol)", "number"),
        ("logp", "Lipophilie logP", "signed_number"),
        ("logd", "Lipophilie logD", "signed_number"),
        ("tpsa", "Polarite TPSA (A^2)", "number"),
        ("hbd", "Donneurs H (HBD)", "count"),
        ("hba", "Accepteurs H (HBA)", "count"),
        ("amide_linker", "Motif amide", "yes_no"),
        ("ester_linker", "Motif ester", "yes_no"),
        ("carbonyl_group", "Groupe carbonyle", "yes_no"),
        ("amine_functionality", "Fonction amine", "text"),
        ("halogenation", "Halogenation", "text"),
        ("fluorination", "Fluoration", "yes_no"),
        ("steric_profile", "Encombrement sterique", "level"),
        ("rotatable_bonds", "Liaisons rotatives", "count"),
        ("ring_constraints", "Contraintes des cycles", "text"),
        ("molecular_shape", "Forme moleculaire (3D)", "text"),
        ("elongation", "Elongation", "text"),
        ("branching", "Ramification", "text"),
        ("chain_length", "Longueur de chaine (C)", "count"),
        ("chain_branching", "Ramification de chaine", "text"),
        ("terminal_group", "Groupe terminal", "text"),
        ("ring_substitution", "Substitution des cycles", "text"),
        ("regiochemistry", "Regiochimie", "text"),
        ("stereochemistry", "Stereochimie", "text"),
        ("conformational_constraint", "Contrainte conformationnelle", "text"),
        ("electronic_profile", "Profil electronique", "text"),
    ]

    PHARMACOLOGY_FIELDS = [
        ("ki_cb1", "Ki CB1 (nM)", "number"),
        ("ki_cb2", "Ki CB2 (nM)", "number"),
        ("ec50", "EC50 (nM)", "number"),
        ("emax", "Emax (%)", "signed_number"),
        ("efficacy_class", "Classe d'activite", "efficacy"),
        ("binding_mode", "Mode de liaison", "text"),
        ("functional_activity", "Activite fonctionnelle", "text"),
        ("assay_system", "Essai / systeme utilise", "text"),
        ("biased_signaling", "Signalisation biaisee observee", "yes_no"),
        ("selectivity", "Selectivite / autres cibles", "text"),
    ]

    PK_FIELDS = [
        ("half_life", "Demi-vie (h)", "number"),
        ("bbb_access_min", "Acces SNC/BHE minimum (%)", "percentage"),
        ("bbb_access_max", "Acces SNC/BHE maximum (%)", "percentage"),
        ("bbb_method", "Mesure BHE / methode", "text"),
        ("membrane_partitioning", "Repartition membranaire", "text"),
        ("protein_binding", "Liaison proteique (%)", "percentage"),
        ("free_fraction", "Fraction libre (%)", "percentage"),
        ("distribution", "Distribution observee", "text"),
        ("metabolic_stability", "Stabilite metabolique", "level"),
        ("active_metabolites", "Metabolites actifs", "text"),
        ("pk_notes", "Notes pharmacocinetiques", "text"),
    ]

    FIELD_GROUPS = {
        "Structure": STRUCTURE_FIELDS,
        "Pharmacologie": PHARMACOLOGY_FIELDS,
        "Pharmacocinetique": PK_FIELDS,
    }

    def __init__(self, master, knowledge):
        super().__init__(master)
        self.knowledge = knowledge
        self.title("Mode simplifie - fiches moleculaires")
        self.geometry("1180x820")
        self.minsize(980, 680)
        self.variables = {}
        self.field_kinds = {}
        self.records = []
        self._create_interface()

    def _create_interface(self):
        heading = ttk.Label(
            self,
            text="Fiches de comparaison experimentale",
            font=("Segoe UI", 16, "bold")
        )
        heading.pack(anchor="w", padx=16, pady=(14, 2))

        ttk.Label(
            self,
            text=(
                "Saisissez des donnees mesurees ou documentees. "
                "Le logiciel ne predit pas l'activite, ne classe pas les molecules "
                "et ne recommande pas de modifications."
            ),
            wraplength=1100
        ).pack(anchor="w", padx=16, pady=(0, 10))

        self.notebook = ttk.Notebook(self)
        self.notebook.pack(fill="both", expand=True, padx=12, pady=6)

        for tab_name, fields in self.FIELD_GROUPS.items():
            tab = ttk.Frame(self.notebook, padding=12)
            self.notebook.add(tab, text=tab_name)
            self._create_field_grid(tab, fields)

        comparison_tab = ttk.Frame(self.notebook, padding=12)
        self.notebook.add(comparison_tab, text="Comparaison")
        self._create_comparison_tab(comparison_tab)

        actions = ttk.Frame(self, padding=(12, 6, 12, 12))
        actions.pack(fill="x")
        ttk.Button(
            actions,
            text="Ajouter a la comparaison",
            command=self._add_record
        ).pack(side="right")

    def _create_field_grid(self, parent, fields):
        panels = []
        for column in range(2):
            panel = ttk.Frame(parent)
            panel.grid(row=0, column=column, sticky="new", padx=(0, 12))
            panel.columnconfigure(1, weight=1)
            parent.columnconfigure(column, weight=1)
            panels.append(panel)

        for index, (key, label, kind) in enumerate(fields):
            panel = panels[index % 2]
            row = index // 2
            variable = tk.StringVar()
            self.variables[key] = variable
            self.field_kinds[key] = kind
            ttk.Label(panel, text=label, width=31).grid(
                row=row, column=0, sticky="w", padx=(0, 8), pady=5
            )

            if kind == "text" or kind in ("number", "signed_number", "percentage", "count"):
                widget = ttk.Entry(panel, textvariable=variable, width=28)
            else:
                widget = ttk.Combobox(
                    panel,
                    textvariable=variable,
                    values=self._choice_values(kind),
                    state="readonly",
                    width=25
                )
                variable.set("Non renseigne")
            widget.grid(row=row, column=1, sticky="ew", pady=5)

    def _choice_values(self, kind):
        if kind == "scaffold":
            examples = self.knowledge["structural_features"]["core_scaffold"]["examples"]
            return ["Non renseigne"] + examples
        if kind == "efficacy":
            values = self.knowledge["metrics"]["efficacy_class"]["values"]
            return ["Non renseigne"] + [value.replace("_", " ") for value in values]
        if kind == "yes_no":
            return ["Non renseigne", "Oui", "Non"]
        if kind == "level":
            return ["Non renseigne", "Faible", "Moyen", "Eleve", "Mesure numerique dans les notes"]
        return []

    def _create_comparison_tab(self, parent):
        columns = ("name", "ki_cb1", "ki_cb2", "half_life", "logp", "bbb")
        self.comparison = ttk.Treeview(parent, columns=columns, show="headings", height=8)
        headings = {
            "name": ("Molecule", 190),
            "ki_cb1": ("Ki CB1 (nM)", 115),
            "ki_cb2": ("Ki CB2 (nM)", 115),
            "half_life": ("Demi-vie (h)", 115),
            "logp": ("logP", 90),
            "bbb": ("Acces SNC/BHE (%)", 170),
        }
        for column, (label, width) in headings.items():
            self.comparison.heading(column, text=label)
            self.comparison.column(column, width=width, anchor="w")
        self.comparison.bind("<<TreeviewSelect>>", self._show_selected_record)
        self.comparison.pack(fill="x", pady=(0, 10))

        ttk.Label(
            parent,
            text="Selectionnez une ligne pour consulter tous les parametres de la fiche."
        ).pack(anchor="w", pady=(0, 8))

        report_frame = ttk.Frame(parent)
        report_frame.pack(fill="both", expand=True)
        self.report = tk.Text(
            report_frame,
            wrap="word",
            height=14,
            font=("Segoe UI", 10),
            padx=10,
            pady=10
        )
        scrollbar = ttk.Scrollbar(report_frame, command=self.report.yview)
        self.report.configure(yscrollcommand=scrollbar.set)
        self.report.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")
        self.report.insert(
            "1.0",
            "Ajoutez une fiche pour comparer les valeurs saisies. "
            "Les mesures restent dans cette fenetre et ne sont pas enregistrees sur disque."
        )
        self.report.configure(state="disabled")

    def _validated_values(self):
        values = {key: variable.get().strip() for key, variable in self.variables.items()}
        if not values["name"]:
            raise ValueError("Saisissez un nom ou un identifiant pour cette molecule.")

        parsed = {}
        for key, value in values.items():
            kind = self.field_kinds[key]
            if not value or value == "Non renseigne" or kind not in (
                "number", "signed_number", "percentage", "count"
            ):
                continue
            try:
                number = float(value.replace(",", "."))
            except ValueError as error:
                raise ValueError(f"La valeur '{value}' n'est pas numerique.") from error
            if not math.isfinite(number):
                raise ValueError("Les valeurs numeriques doivent etre finies.")
            if kind != "signed_number" and number < 0:
                raise ValueError("Les valeurs de comptage et de concentration doivent etre positives.")
            if kind == "count" and not number.is_integer():
                raise ValueError("Les nombres d'atomes et de cycles doivent etre entiers.")
            if key in ("ki_cb1", "ki_cb2", "ec50", "half_life") and number <= 0:
                raise ValueError("Ki, EC50 et demi-vie doivent etre superieurs a zero.")
            if kind == "percentage" and number > 100:
                raise ValueError("Une valeur en pourcentage doit etre comprise entre 0 et 100.")
            parsed[key] = number

        minimum = parsed.get("bbb_access_min")
        maximum = parsed.get("bbb_access_max")
        if minimum is not None and maximum is not None and minimum > maximum:
            raise ValueError("La borne minimale d'acces BHE depasse la borne maximale.")

        return values, parsed

    def _add_record(self):
        try:
            values, parsed = self._validated_values()
            descriptors = describe_smiles(values["smiles"]) if values["smiles"] else None
        except ImportError:
            messagebox.showerror(
                "Dependance manquante",
                "Installez les dependances avec : python -m pip install -r requirements.txt",
                parent=self
            )
            return
        except ValueError as error:
            messagebox.showerror("Fiche invalide", str(error), parent=self)
            return

        record = {**values, "_parsed": parsed, "_descriptors": descriptors}
        self.records.append(record)
        self.comparison.insert(
            "",
            "end",
            iid=str(len(self.records)),
            values=(
                values["name"],
                values["ki_cb1"] or "-",
                values["ki_cb2"] or "-",
                values["half_life"] or "-",
                values["logp"] or "-",
                self._format_bbb(values)
            )
        )
        self._show_record(record)
        self.notebook.select(3)

    @staticmethod
    def _format_bbb(values):
        minimum = values["bbb_access_min"]
        maximum = values["bbb_access_max"]
        if minimum and maximum:
            return f"{minimum}-{maximum}"
        if minimum:
            return f">= {minimum}"
        if maximum:
            return f"<= {maximum}"
        return "-"

    def _show_record(self, record):
        sections = [
            ("Structure", self.STRUCTURE_FIELDS),
            ("Pharmacologie", self.PHARMACOLOGY_FIELDS),
            ("Pharmacocinetique", self.PK_FIELDS),
        ]
        lines = [
            f"Fiche: {record['name']}",
            "Donnees saisies: experimentales ou documentees, non predites.",
            "Aucun score, classement ou conseil pharmacologique n'est produit.",
        ]

        for title, fields in sections:
            lines.extend(["", title])
            for key, label, _kind in fields:
                value = record[key]
                if value and value != "Non renseigne":
                    lines.append(f"- {label}: {value}")

        if record["_descriptors"]:
            descriptors = record["_descriptors"]
            lines.extend([
                "",
                "Descripteurs calcules depuis le SMILES (RDKit)",
                f"- Formule: {descriptors['formula']}",
                f"- Masse exacte: {descriptors['molecular_weight']:.3f} g/mol",
                f"- cLogP calcule: {descriptors['clogp']:.2f}",
                f"- TPSA: {descriptors['tpsa']:.2f} A^2",
                f"- HBD / HBA: {descriptors['hbd']} / {descriptors['hba']}",
                f"- Cycles / cycles aromatiques: {descriptors['ring_count']} / {descriptors['aromatic_ring_count']}",
                f"- Liaisons rotatives: {descriptors['rotatable_bonds']}",
            ])

        lines.extend([
            "",
            "Contexte et limites du dataset",
            "- Ki, EC50 et Emax mesurent des proprietes differentes; Ki ne permet pas de deduire EC50 ou Emax.",
            "- Les observations SAR dependent des series et des essais; elles ne sont pas des regles universelles.",
            "- La permeabilite BHE ne constitue pas une mesure d'affinite CB1.",
            "- La demi-vie est une donnee saisie par l'utilisateur, pas une propriete calculee par le dataset.",
        ])

        self.report.configure(state="normal")
        self.report.delete("1.0", tk.END)
        self.report.insert("1.0", "\n".join(lines))
        self.report.configure(state="disabled")

    def _show_selected_record(self, event=None):
        selection = self.comparison.selection()
        if selection:
            record_index = int(selection[0]) - 1
            self._show_record(self.records[record_index])