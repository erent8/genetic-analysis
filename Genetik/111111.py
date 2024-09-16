import tkinter as tk
from tkinter import ttk, scrolledtext
from Bio import Entrez
from Bio import SeqIO

# NCBI Entrez e-posta adresini tanımla
Entrez.email = "themimar080@gmail.com"

def fetch_pdyn_gene(accession_id):
    try:
        # GenBank verilerini NCBI'dan çek
        handle = Entrez.efetch(db="nucleotide", id=accession_id, rettype="gb", retmode="text")
        record = SeqIO.read(handle, "genbank")
        handle.close()

        # Organizma bilgisi ve gen uzunluğu
        output = []
        output.append("===== GENEL BİLGİLER =====")
        output.append(f"Organizma: {record.annotations['organism']}")
        output.append(f"Dizi Uzunluğu: {len(record.seq)} baz çifti")
        
        # Genel açıklamalar ve meta veriler
        output.append("\n===== ANNOTATIONS (AÇIKLAMALAR) =====")
        for key, value in record.annotations.items():
            output.append(f"{key.capitalize()}: {value}")

        # Tüm feature'lar
        output.append("\n===== GENİN ÖZELLİKLERİ (FEATURES) =====")
        for feature in record.features:
            output.append(f"Özellik Türü: {feature.type}")
            output.append(f"Konum: {feature.location}")
            if feature.qualifiers:
                for qual_key, qual_value in feature.qualifiers.items():
                    output.append(f"{qual_key.capitalize()}: {qual_value}")
            output.append("-" * 40)  # Feature'lar arasına çizgi ekliyoruz

        # CDS (Protein) bilgisi
        for feature in record.features:
            if feature.type == "CDS":
                output.append("\n===== PROTEİN BİLGİLERİ =====")
                output.append(f"Ürün: {feature.qualifiers.get('product', ['Bilinmiyor'])[0]}")
                output.append(f"Protein ID: {feature.qualifiers.get('protein_id', ['Bilinmiyor'])[0]}")
                output.append(f"Amino Asit Dizisi (İlk 10): {feature.qualifiers.get('translation', ['Bilinmiyor'])[0][:10]}...")
                break  # Sadece ilk CDS'yi gösteriyoruz

        # DNA dizisinin ilk 100 baz çifti
        output.append("\n===== DNA DİZİSİ (İlk 100 Baz Çifti) =====")
        sequence_str = str(record.seq[:100])  # Diziyi string'e çeviriyoruz
        for i in range(0, len(sequence_str), 50):  # Her 50 bazdan sonra satır atla
            output.append(sequence_str[i:i+50])
        
        return "\n".join(output)

    except Exception as e:
        return f"Veri çekilirken hata oluştu: {e}"

# Tkinter GUI oluşturma
def create_gui(data):
    root = tk.Tk()
    root.title("PDYN Gen Bilgileri")
    
    # Pencere boyutlarını ayarla
    root.geometry("800x600")
    
    # ScrolledText widget'ını oluştur
    text = scrolledtext.ScrolledText(root, wrap=tk.WORD, font=("Arial", 12))
    text.pack(fill=tk.BOTH, expand=True)
    
    # Verileri ScrolledText widget'ına ekle
    text.insert(tk.END, data)
    
    # Uygulamayı başlat
    root.mainloop()

# PDYN geni için Homo sapiens erişim numarası
pdyn_accession_id = "NM_024411"

# Verileri çek ve GUI'de göster
gene_data = fetch_pdyn_gene(pdyn_accession_id)
create_gui(gene_data)
