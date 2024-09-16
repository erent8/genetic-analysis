# api_address "a37830e44ae146dfeb799dfc589066216a08"
# profiles_id = "themimar080@gmail.com"

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
        print("\n===== GENEL BİLGİLER =====")
        print(f"Organizma: {record.annotations['organism']}")
        print(f"Dizi Uzunluğu: {len(record.seq)} baz çifti")
        
        # Genel açıklamalar ve meta veriler
        print("\n===== ANNOTATIONS (AÇIKLAMALAR) =====")
        for key, value in record.annotations.items():
            print(f"{key.capitalize()}: {value}")

        # Tüm feature'lar
        print("\n===== GENİN ÖZELLİKLERİ (FEATURES) =====")
        for feature in record.features:
            print(f"Özellik Türü: {feature.type}")
            print(f"Konum: {feature.location}")
            if feature.qualifiers:
                for qual_key, qual_value in feature.qualifiers.items():
                    print(f"{qual_key.capitalize()}: {qual_value}")
            print("-" * 40)  # Feature'lar arasına çizgi ekliyoruz

        # CDS (Protein) bilgisi
        for feature in record.features:
            if feature.type == "CDS":
                print("\n===== PROTEİN BİLGİLERİ =====")
                print(f"Ürün: {feature.qualifiers.get('product', ['Bilinmiyor'])[0]}")
                print(f"Protein ID: {feature.qualifiers.get('protein_id', ['Bilinmiyor'])[0]}")
                print(f"Amino Asit Dizisi (İlk 10): {feature.qualifiers.get('translation', ['Bilinmiyor'])[0][:10]}...")
                break  # Sadece ilk CDS'yi gösteriyoruz

        # DNA dizisinin ilk 100 baz çifti
        print("\n===== DNA DİZİSİ (İlk 100 Baz Çifti) =====")
        sequence_str = str(record.seq[:100])  # Diziyi string'e çeviriyoruz
        for i in range(0, len(sequence_str), 50):  # Her 50 bazdan sonra satır atla
            print(sequence_str[i:i+50])
        
    except Exception as e:
        print(f"Veri çekilirken hata oluştu: {e}")

# PDYN geni için Homo sapiens erişim numarası
pdyn_accession_id = "NM_024411"

# Verileri çek ve terminale yazdır
fetch_pdyn_gene(pdyn_accession_id)
