import tkinter as tk
from tkinter import ttk, messagebox
from sar_engine import build_knowledge_report, describe_smiles, load_knowledge_base
from simple_mode import SimpleMoleculeTester

SCAFFOLDS = {
    "Indole": {
        "nitrogen_count": 1,
        "ring_count": 2,
        "aromatic": True
    },
    "Indazole": {
        "nitrogen_count": 2,
        "ring_count": 2,
        "aromatic": True
    }
}

FUNCTIONAL_GROUPS = {
    "Amides": {
        "Carboxamide primaire": "R-C(=O)-NH2",
        "Carboxamide secondaire": "R-C(=O)-NHR",
        "Carboxamide tertiaire": "R-C(=O)-NR2",
        "Lactame": "Amide cyclique",
        "Thioamide": "R-C(=S)-NR2"
    },

    "Esters": {
        "Ester carboxylique": "R-C(=O)-OR",
        "Lactone": "Ester cyclique",
        "Thioester": "R-C(=O)-SR"
    },

    "Acides carboxyliques": {
        "Acide carboxylique": "R-C(=O)-OH"
    },

    "Amines": {
        "Amine primaire": "R-NH2",
        "Amine secondaire": "R2-NH",
        "Amine tertiaire": "R3-N"
    },

    "Alcools et phénols": {
        "Alcool": "R-OH",
        "Phénol": "Ar-OH"
    },

    "Éthers": {
        "Éther": "R-O-R",
        "Éther cyclique": "O dans un cycle"
    },

    "Composés carbonylés": {
        "Aldéhyde": "R-C(=O)-H",
        "Cétone": "R-C(=O)-R"
    },

    "Fonctions azotées": {
        "Nitrile": "R-C≡N",
        "Nitro": "R-NO2",
        "Hydrazine": "R-NH-NH-R"
    },

    "Fonctions soufrées": {
        "Thiol": "R-SH",
        "Thioéther": "R-S-R",
        "Thioamide": "R-C(=S)-NR2",
        "Thioester": "R-C(=O)-SR"
    },

    "Halogénures": {
        "Fluorure": "R-F",
        "Chlorure": "R-Cl",
        "Bromure": "R-Br",
        "Iodure": "R-I"
    },

    "Liaisons multiples": {
        "Alcène": "C=C",
        "Alcyne": "C≡C"
    }
}

class MolecularAnalyzer:

    def __init__(self, nitrogen_count, ring_count, aromatic):
        self.nitrogen_count = nitrogen_count
        self.ring_count = ring_count
        self.aromatic = aromatic

    def identify_scaffold(self):
        candidates = []

        for name, data in SCAFFOLDS.items():

            if data["nitrogen_count"] != self.nitrogen_count:
                continue

            if data["ring_count"] != self.ring_count:
                continue

            if data["aromatic"] != self.aromatic:
                continue

            candidates.append(name)

        return candidates

    def analyze(self):
        candidates = self.identify_scaffold()

        if not candidates:
            return "Aucun squelette connu ne correspond exactement aux paramètres."

        return "Squelette(s) correspondant(s) :\n\n" + "\n".join(
            f"• {x}" for x in candidates
        )


class Application(tk.Tk):

    def __init__(self):
        super().__init__()

        self.title("Analyseur structural et SAR")
        self.geometry("850x700")
        self.minsize(750, 600)

        self.configure(bg="#202020")
        self.knowledge = load_knowledge_base()

        self.create_interface()

    def create_interface(self):

        title = tk.Label(
            self,
            text="Analyseur structural et SAR",
            font=("Segoe UI", 22, "bold"),
            bg="#202020",
            fg="white"
        )
        title.pack(pady=(20, 5))

        subtitle = tk.Label(
            self,
            text="Classification de structures et groupes fonctionnels",
            font=("Segoe UI", 11),
            bg="#202020",
            fg="#bbbbbb"
        )
        subtitle.pack(pady=(0, 20))

        ttk.Button(
            self,
            text="Ouvrir le mode simplifie",
            command=self.open_simple_mode
        ).pack(pady=(0, 8))

        main = tk.Frame(
            self,
            bg="#292929"
        )
        main.pack(
            fill="both",
            expand=True,
            padx=25,
            pady=10
        )

        self.create_structure_section(main)
        self.create_functional_group_section(main)
        self.create_result_section(main)

    def open_simple_mode(self):
        SimpleMoleculeTester(self, self.knowledge)

    def create_structure_section(self, parent):

        frame = tk.LabelFrame(
            parent,
            text="Paramètres structuraux",
            font=("Segoe UI", 11, "bold"),
            bg="#292929",
            fg="white",
            padx=15,
            pady=15
        )

        frame.pack(
            fill="x",
            padx=15,
            pady=15
        )

        tk.Label(
            frame,
            text="Azotes cycliques :",
            bg="#292929",
            fg="white"
        ).grid(row=0, column=0, sticky="w", pady=5)

        self.nitrogen_var = tk.StringVar(value="1")

        tk.Label(
            frame,
            textvariable=self.nitrogen_var,
            bg="#292929",
            fg="white"
        ).grid(row=0, column=1, padx=15)

        tk.Label(
            frame,
            text="Nombre de cycles :",
            bg="#292929",
            fg="white"
        ).grid(row=1, column=0, sticky="w", pady=5)

        self.ring_var = tk.StringVar(value="2")

        tk.Label(
            frame,
            textvariable=self.ring_var,
            bg="#292929",
            fg="white"
        ).grid(row=1, column=1, padx=15)

        tk.Label(
            frame,
            text="Aromaticité :",
            bg="#292929",
            fg="white"
        ).grid(row=2, column=0, sticky="w", pady=5)

        self.aromatic_var = tk.StringVar(value="Oui")

        tk.Label(
            frame,
            textvariable=self.aromatic_var,
            bg="#292929",
            fg="white"
        ).grid(row=2, column=1, sticky="w", padx=15)

        tk.Label(
            frame,
            text="Molécule (SMILES) :",
            bg="#292929",
            fg="white"
        ).grid(row=3, column=0, sticky="w", pady=5)

        self.smiles_var = tk.StringVar(value="c1ccc2[nH]ccc2c1")

        tk.Entry(
            frame,
            textvariable=self.smiles_var,
            width=35
        ).grid(row=3, column=1, padx=15, sticky="ew")

        ttk.Button(
            frame,
            text="Visualiser en 3D",
            command=self.visualize_molecule
        ).grid(row=3, column=2, padx=8)

    def visualize_molecule(self):
        try:
            from rdkit import Chem
            from rdkit.Chem import AllChem
            from matplotlib.backends.backend_tkagg import (
                FigureCanvasTkAgg,
                NavigationToolbar2Tk
            )
            from matplotlib.figure import Figure
        except ImportError:
            messagebox.showerror(
                "Dépendances manquantes",
                "Installez les dépendances avec :\n"
                "python -m pip install -r requirements.txt"
            )
            return

        smiles = self.smiles_var.get().strip()
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            messagebox.showerror(
                "SMILES invalide",
                "Vérifiez la formule SMILES saisie."
            )
            return

        molecule = Chem.AddHs(molecule)
        if AllChem.EmbedMolecule(molecule, randomSeed=42) != 0:
            messagebox.showerror(
                "Géométrie indisponible",
                "La conformation 3D n'a pas pu être calculée."
            )
            return

        if AllChem.MMFFHasAllMoleculeParams(molecule):
            AllChem.MMFFOptimizeMolecule(molecule)
        else:
            AllChem.UFFOptimizeMolecule(molecule)

        positions = molecule.GetConformer().GetPositions()
        atom_indices = [
            atom.GetIdx() for atom in molecule.GetAtoms()
            if atom.GetSymbol() != "H"
        ]
        colors = {
            "C": "#555b64", "N": "#3978d4", "O": "#d94a45",
            "S": "#d6a91e", "P": "#e27d32", "F": "#55a65a",
            "Cl": "#55a65a", "Br": "#8e5bb7", "I": "#75418c"
        }

        window = tk.Toplevel(self)
        window.title("Visualisation moléculaire 3D")
        window.geometry("760x620")

        figure = Figure(figsize=(7, 5), dpi=100, facecolor="#f4f5f6")
        axes = figure.add_subplot(111, projection="3d")
        axes.set_facecolor("#f4f5f6")

        for bond in molecule.GetBonds():
            start = bond.GetBeginAtomIdx()
            end = bond.GetEndAtomIdx()
            if start not in atom_indices or end not in atom_indices:
                continue
            axes.plot(
                positions[[start, end], 0],
                positions[[start, end], 1],
                positions[[start, end], 2],
                color="#737b84",
                linewidth=2.2
            )

        for index in atom_indices:
            atom = molecule.GetAtomWithIdx(index)
            x, y, z = positions[index]
            axes.scatter(
                x, y, z,
                color=colors.get(atom.GetSymbol(), "#b34b8c"),
                s=170,
                edgecolors="white",
                linewidths=0.8
            )
            axes.text(
                x, y, z, atom.GetSymbol(), color="white",
                ha="center", va="center", fontsize=8, fontweight="bold"
            )

        center = positions[atom_indices].mean(axis=0)
        radius = max(
            (positions[atom_indices].max(axis=0) -
             positions[atom_indices].min(axis=0)).max() / 2,
            1.5
        )
        axes.set_xlim(center[0] - radius, center[0] + radius)
        axes.set_ylim(center[1] - radius, center[1] + radius)
        axes.set_zlim(center[2] - radius, center[2] + radius)
        axes.set_box_aspect((1, 1, 1))
        axes.set_axis_off()
        axes.set_title(smiles, fontsize=10, pad=12)

        canvas = FigureCanvasTkAgg(figure, master=window)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        toolbar = NavigationToolbar2Tk(canvas, window, pack_toolbar=False)
        toolbar.update()
        toolbar.pack(side="bottom", fill="x")

    def create_functional_group_section(self, parent):

        frame = tk.LabelFrame(
            parent,
            text="Groupe fonctionnel",
            font=("Segoe UI", 11, "bold"),
            bg="#292929",
            fg="white",
            padx=15,
            pady=15
        )

        frame.pack(
            fill="x",
            padx=15,
            pady=10
        )

        tk.Label(
            frame,
            text="Famille :",
            bg="#292929",
            fg="white"
        ).grid(row=0, column=0, sticky="w", pady=5)

        self.family_var = tk.StringVar()

        families = list(FUNCTIONAL_GROUPS.keys())

        self.family_combo = ttk.Combobox(
            frame,
            textvariable=self.family_var,
            values=families,
            state="readonly",
            width=30
        )

        self.family_combo.grid(
            row=0,
            column=1,
            padx=15
        )

        self.family_combo.bind(
            "<<ComboboxSelected>>",
            self.update_subtypes
        )

        tk.Label(
            frame,
            text="Sous-type :",
            bg="#292929",
            fg="white"
        ).grid(row=1, column=0, sticky="w", pady=5)

        self.subtype_var = tk.StringVar()

        self.subtype_combo = ttk.Combobox(
            frame,
            textvariable=self.subtype_var,
            state="readonly",
            width=30
        )

        self.subtype_combo.grid(
            row=1,
            column=1,
            padx=15
        )

        self.pattern_label = tk.Label(
            frame,
            text="",
            font=("Consolas", 11),
            bg="#292929",
            fg="#66ccff"
        )

        self.pattern_label.grid(
            row=2,
            column=0,
            columnspan=2,
            pady=15
        )

        self.subtype_combo.bind(
            "<<ComboboxSelected>>",
            self.update_pattern
        )

    def update_subtypes(self, event=None):

        family = self.family_var.get()

        if not family:
            return

        subtypes = list(
            FUNCTIONAL_GROUPS[family].keys()
        )

        self.subtype_combo["values"] = subtypes

        if subtypes:
            self.subtype_var.set(subtypes[0])
            self.update_pattern()

    def update_pattern(self, event=None):

        family = self.family_var.get()
        subtype = self.subtype_var.get()

        if not family or not subtype:
            return

        pattern = FUNCTIONAL_GROUPS[family][subtype]

        self.pattern_label.config(
            text=f"Structure : {pattern}"
        )

    def create_result_section(self, parent):

        frame = tk.LabelFrame(
            parent,
            text="Analyse",
            font=("Segoe UI", 11, "bold"),
            bg="#292929",
            fg="white",
            padx=15,
            pady=15
        )

        frame.pack(
            fill="both",
            expand=True,
            padx=15,
            pady=10
        )

        self.result_text = tk.Text(
            frame,
            bg="#181818",
            fg="white",
            insertbackground="white",
            font=("Consolas", 11),
            wrap="word"
        )

        scrollbar = ttk.Scrollbar(frame, command=self.result_text.yview)
        self.result_text.configure(yscrollcommand=scrollbar.set)
        self.result_text.pack(side="left", fill="both", expand=True)
        scrollbar.pack(side="right", fill="y")

        button = tk.Button(
            frame,
            text="Analyser",
            command=self.analyze,
            font=("Segoe UI", 11, "bold"),
            padx=20,
            pady=8
        )

        button.pack(
            pady=15
        )

    def analyze(self):
        try:
            descriptors = describe_smiles(self.smiles_var.get())
        except ImportError:
            messagebox.showerror(
                "Dependances manquantes",
                "Installez les dependances avec :\n"
                "python -m pip install -r requirements.txt"
            )
            return
        except ValueError as error:
            messagebox.showerror(
                "SMILES invalide",
                str(error)
            )
            return

        self.nitrogen_var.set(str(descriptors["ring_nitrogen_count"]))
        self.ring_var.set(str(descriptors["ring_count"]))
        self.aromatic_var.set("Oui" if descriptors["aromatic"] else "Non")

        analyzer = MolecularAnalyzer(
            descriptors["ring_nitrogen_count"],
            descriptors["ring_count"],
            descriptors["aromatic"]
        )

        result = (
            "Correspondances de squelette possibles (comptes structurels uniquement) :\n\n"
            + analyzer.analyze()
        )

        family = self.family_var.get()
        subtype = self.subtype_var.get()

        if family and subtype:

            pattern = FUNCTIONAL_GROUPS[family][subtype]

            result += (
                "\n\n"
                "Groupe fonctionnel sélectionné :\n"
                f"• Famille : {family}\n"
                f"• Sous-type : {subtype}\n"
                f"• Motif : {pattern}"
            )

        result += "\n\n" + build_knowledge_report(descriptors, self.knowledge)

        self.result_text.delete(
            "1.0",
            tk.END
        )

        self.result_text.insert(
            tk.END,
            result
        )


if __name__ == "__main__":
    app = Application()
    app.mainloop()