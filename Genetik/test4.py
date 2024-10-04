import tkinter as tk
from tkinter import ttk
from Bio import Entrez
from Bio import SeqIO

# NCBI Entrez
Entrez.email = "YOUR_EMAİL"

def fetch_pdyn_gene(accession_id):
    try:
        # GenBank verilerini NCBI'dan çek
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Organizma bilgisi ve genel verileri döndür
        gene_data = []
        gene_data.append(("Organism", record.annotations['organism']))
        gene_data.append(("Sequence Length", f"{len(record.seq)} bp"))

        # Annotations verilerini ekle
        gene_data.append(("Molecule Type", record.annotations.get('molecule_type', 'N/A')))
        gene_data.append(("Topology", record.annotations.get('topology', 'N/A')))
        gene_data.append(("Date", record.annotations.get('date', 'N/A')))

        # Genin özelliklerini ekle
        for feature in record.features:
            gene_data.append((f"Feature: {feature.type}", f"Location: {feature.location}"))
            for key, value in feature.qualifiers.items():
                gene_data.append((key.capitalize(), ", ".join(value)))

        # CDS bilgilerini ekle (Protein Bilgisi)
        for feature in record.features:
            if feature.type == "CDS":
                gene_data.append(("Protein Product", feature.qualifiers.get('product', ['N/A'])[0]))
                gene_data.append(("Protein ID", feature.qualifiers.get('protein_id', ['N/A'])[0]))
                gene_data.append(("Translation", feature.qualifiers.get('translation', ['N/A'])[0][:10] + "..."))
                break

        return gene_data

    except Exception as e:
        print(f"Error fetching data: {e}")
        return []

# Tkinter GUI oluşturma
def create_gui(data):
    root = tk.Tk()
    root.title("PDYN Gen Bilgileri")
    
    # Pencere boyutlarını ayarla
    root.geometry("800x600")
    
    # Treeview tablosunu oluştur
    tree = ttk.Treeview(root, columns=("Feature", "Value"), show='headings', height=25)
    tree.heading("Feature", text="Özellik")
    tree.heading("Value", text="Değer")
    
    # Sütun genişliklerini ayarla
    tree.column("Feature", width=250, anchor=tk.W)
    tree.column("Value", width=550, anchor=tk.W)

    # Tabloya çizgiler ekle
    style = ttk.Style()
    style.configure("Treeview",
                    background="white",
                    foreground="black",
                    rowheight=30,
                    fieldbackground="lightgray")
    style.configure("Treeview.Heading",
                    background="gray",
                    foreground="white",
                    font=("Arial", 10, "bold"))
    style.map("Treeview.Heading",
              background=[('selected', 'darkgray')])
    
    # Verileri Treeview'e ekle
    for item in data:
        tree.insert("", "end", values=item)
    
    # Treeview'i pencereye yerleştir
    tree.pack(fill=tk.BOTH, expand=True)
    
    # Kaydırma çubukları ekle
    vsb = ttk.Scrollbar(root, orient="vertical", command=tree.yview)
    vsb.pack(side='right', fill='y')
    tree.configure(yscrollcommand=vsb.set)

    hsb = ttk.Scrollbar(root, orient="horizontal", command=tree.xview)
    hsb.pack(side='bottom', fill='x')
    tree.configure(xscrollcommand=hsb.set)

    # Uygulamayı başlat
    root.mainloop()

# PDYN geni için Homo sapiens erişim numarası
pdyn_accession_id = "NM_024411"

# GenBank verilerini çek ve GUI'de göster
gene_data = fetch_pdyn_gene(pdyn_accession_id)
create_gui(gene_data)
