from Bio import SeqIO, AlignIO
import pandas as pd
from collections import Counter
import os


records = []
for record in SeqIO.parse("./analysis/MAfile.fa", "fasta"):
    if 'X' in str(record.seq)[1068:1184] or 'B' in str(record.seq)[1068:1184] or 'J' in str(record.seq)[1068:1184] or 'Z' in str(record.seq)[1068:1184]: # filter out missing information in region of interest
        continue
    else:
        records.append(record)

SeqIO.write(records,'./analysis/MAfile_filt.fa','fasta')



records_dengue_1 = [] ; records_dengue_2 = [] ; records_dengue_3 = [] ; records_dengue_4 = [] ; records_Japanese_encephalitis = [] 
records_Kyasanur_Forest_disease = [] ; records_Langat = [] ; records_Louping_ill = [] ; records_Murray_Valley_encephalitis = [] ; records_Omsk_hemorrhagic_fever = [] 
records_Powassan = [] ; records_Saint_Louis_encephalitis = [] ; records_Tick_borne_encephalitis = [] ; records_Usutu = [] ; records_West_Nile = [] ; records_Yellow_fever = [] ; records_Zika = []
for record in SeqIO.parse("./analysis/MAfile_filt.fa", "fasta"):
    if any([x in record.id for x in ['dengue_virus_type_1','dengue_virus_type_I','Dengue_virus_1']]):
        records_dengue_1.append(record)
    elif any([x in record.id for x in ['Dengue_virus_type_2','dengue_virus_type_2','Dengue_virus_2']]):
        records_dengue_2.append(record)
    elif any([x in record.id for x in ['Dengue_virus_type_3','dengue_virus_type_3','Dengue_virus_3']]):
        records_dengue_3.append(record)
    elif any([x in record.id for x in ['Dengue_virus_type_4','dengue_virus_type_4','Dengue_virus_4']]):
        records_dengue_4.append(record)
    elif 'Japanese_encephalitis_virus' in record.id:
        records_Japanese_encephalitis.append(record)
    elif 'Kyasanur_Forest_disease_virus' in record.id:
        records_Kyasanur_Forest_disease.append(record)
    elif 'Langat_virus' in record.id:
        records_Langat.append(record)
    elif 'Louping_ill_virus' in record.id:
        records_Louping_ill.append(record)
    elif 'Murray_Valley_encephalitis_virus' in record.id:
        records_Murray_Valley_encephalitis.append(record)
    elif 'Omsk_hemorrhagic_fever_virus' in record.id:
        records_Omsk_hemorrhagic_fever.append(record)
    elif 'Saint_Louis_encephalitis_virus' in record.id:
        records_Saint_Louis_encephalitis.append(record)
    elif 'Tick-borne_encephalitis_virus' in record.id:
        records_Tick_borne_encephalitis.append(record)
    elif 'Usutu_virus' in record.id:
        records_Usutu.append(record)
    elif 'West_Nile_virus' in record.id:
        records_West_Nile.append(record)
    elif 'Yellow_fever_virus' in record.id:
        records_Yellow_fever.append(record)
    elif any([x in record.id for x in ['Powassan_virus','Deer_tick_virus']]):
        records_Powassan.append(record)
    elif 'Zika_virus' in record.id:
        records_Zika.append(record)


SeqIO.write(records_dengue_1,'./analysis/separated_strains/Dengue_virus_type_1.fa','fasta')
SeqIO.write(records_dengue_2,'./analysis/separated_strains/Dengue_virus_type_2.fa','fasta')
SeqIO.write(records_dengue_3,'./analysis/separated_strains/Dengue_virus_type_3.fa','fasta')
SeqIO.write(records_dengue_4,'./analysis/separated_strains/Dengue_virus_type_4.fa','fasta')
SeqIO.write(records_Japanese_encephalitis,'./analysis/separated_strains/Japanese_encephalitis.fa','fasta')
SeqIO.write(records_Kyasanur_Forest_disease,'./analysis/separated_strains/Kyasanur_Forest_disease.fa','fasta')
SeqIO.write(records_Langat,'./analysis/separated_strains/Langat.fa','fasta')
SeqIO.write(records_Louping_ill,'./analysis/separated_strains/Louping_ill.fa','fasta')
SeqIO.write(records_Murray_Valley_encephalitis,'./analysis/separated_strains/Murray_Valley_encephalitis.fa','fasta')
SeqIO.write(records_Omsk_hemorrhagic_fever,'./analysis/separated_strains/Omsk_hemorrhagic_fever.fa','fasta')
SeqIO.write(records_Saint_Louis_encephalitis,'./analysis/separated_strains/Saint_Louis_encephalitis.fa','fasta')
SeqIO.write(records_Tick_borne_encephalitis,'./analysis/separated_strains/Tick_borne_encephalitis.fa','fasta')
SeqIO.write(records_Usutu,'./analysis/separated_strains/Usutu.fa','fasta')
SeqIO.write(records_West_Nile,'./analysis/separated_strains/West_Nile.fa','fasta')
SeqIO.write(records_Yellow_fever,'./analysis/separated_strains/Yellow_fever.fa','fasta')
SeqIO.write(records_Powassan,'./analysis/separated_strains/Powassan.fa','fasta')
SeqIO.write(records_Zika,'./analysis/separated_strains/Zika.fa','fasta')

counter = 0
files = [filename for filename in os.listdir('./analysis/separated_strains')]
for file in files:
    counter += 1
    print('calculating aa frequencies for ' + file.replace('.fa', '') + ' [' + str(counter) +' of ' + str(len(files))+']')
    alignment = AlignIO.read("./analysis/separated_strains/" + file, "fasta")
    aaFreqs = {}
    for a in range(len(alignment[0].seq)):
        aaFreqs[a] = Counter(alignment[:, a])
    df = pd.DataFrame(aaFreqs) ; df_T = df.T
    df_T.to_csv('./analysis/consensus/' + file.replace('.fa','.csv'))
