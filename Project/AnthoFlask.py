from flask import Flask, render_template, request
from Bio import Entrez
from Bio import Medline
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize
import re

app = Flask(__name__)


def get_number_of_hits(mail, query):
    Entrez.email = mail
    handle = Entrez.egquery(term=query)
    record = Entrez.read(handle)
    for row in record["eGQueryResult"]:
        if row["DbName"] == "pubmed":
            maxval = (row["Count"])

            return maxval


def get_records(mail, query,
                hits):  # methode voor het ophalen van abstracts uit pubmed die horen bij de query. input=mail, query, hits. output=record (dit zijn de abstracts)
    Entrez.email = mail
    handle = Entrez.esearch(db="pubmed", term=query, retmax=hits)
    record = Entrez.read(handle)
    idlist = record["IdList"]
    handle = Entrez.efetch(db="pubmed", id=idlist, rettype="medline", retmode="text")
    records = Medline.parse(handle)

    return records


def process_records(
        records):  # methode voor het filteren van pubmed abstracts. input(records) output(info_dict,no_abstract)
    stop_words = set(stopwords.words(
        "english"))  # lijst met meestgebruikte Engelse woorden (en vaak vorkomende wetenschappelijke termen) die uit het abstract gefilterd zullen worden
    most_used_words = ["a", "about", "all", "also", "and", "as", "at", "be", "because", "but", "by", "can", "come",
                       "could", "day", "do", "even", "find", "first", "for", "from", "get", "give", "go", "have", "he",
                       "her", "here", "him", "his", "how", "I", "if", "in", "into", "it", "its", "just", "know", "like",
                       "look", "make", "man", "many", "me", "more", "my", "new", "no", "not", "now", "of", "on", "one",
                       "only", "or", "other", "our", "out", "people", "say", "see", "she", "so", "some", "take", "tell",
                       "than", "that", "the", "their", "them", "then", "there", "these", "they", "thing", "think",
                       "this", "those", "time", "to", "two", "up", "use", "very", "want", "way", "we", "well", "what",
                       "when", "which", "who", "will", "with", "would", "year", "you", "your", ",", ".", ":", ";",
                       "METHODS", "MATERIAL", "SCOPE", "CHEMICAL", "COMPOUNDS", "CONCLUSION", "CID", "AND", "RESULTS",
                       "BACKGROUND", "CONCLUSION", "CONCLUSIONS", "MAIN", ")", "(", "[", "]", "DNA", "RNA", "mRNA",
                       "RNAs", "DISCUSSION", "PubChem", "MATERIALS", "CONTEXT", "MESSAGE", "REVIEW", "KEY", "CO2",
                       "NaCl", "mM", "RNA-Seq", "PCR", "SNP", "QTLs"]
    info_dict = {}  # een dictionary. de key van deze dictionary is een accessiecode met als value een gefilterd abstract
    count = 0
    no_abstract = []  # lijst van accessiecode die geen abstract bevatten

    for record in records:
        unfiltered_ab = []
        filtered_ab = []
        count += 1
        try:
            unfiltered_ab = word_tokenize(record["AB"])
        except KeyError:
            no_abstract.append("PMID:" + str(record["PMID"]))

        for word in unfiltered_ab:
            if word not in stop_words and word not in most_used_words:  # woorden uit most_used_words woden uit het abstract weg gefilterd
                filtered_ab.append(word)
        if filtered_ab:
            info_dict[record["PMID"]] = filtered_ab

    return [info_dict, no_abstract]


def search_genes(
        info_dict):  # Methode voor het vinden van gen namen middels een regular expression input=info_dict output=gene_dict, no_genes
    regex = re.compile(r"[A-z]+\d*[A-Z]")
    no_genes = []  # lijst met accessiecodes van abstracts die volgens de regular expression geen genen bevatten
    gene_dict = {}

    for pmid, abstract in info_dict.items():
        genelist = []

        for word in abstract:
            if regex.match(word):
                if word not in genelist and word.startswith("Pub") is not True:
                    genelist.append(word)
        if genelist:
            gene_dict[pmid] = [genelist, abstract]
        else:
            no_genes.append("PMID:" + pmid)

    return [gene_dict, no_genes]


def search_stress(
        gene_dict):  # Methode voor het vinden van stess omstandigheden middels een regular expression input=gene_dict output=result_dict, no_stress
    no_stress = []  # lijst met accessiecodes van abstracts die volgens de regular expression geen stess omstandigheden bevatten. Dit is de output van deze methode
    result_dict = {}  # dictionary met de volgende opbouw key=accessiecode value=[lijst van gevonden genen in het abstract,lijst van gevonden stresomstandigheden in het abstract]. Dit is de output van deze methode.
    for pmid in gene_dict:
        stresslist = []
        for word in gene_dict[pmid][1]:
            if "stress" in word:
                stresslist.append(gene_dict[pmid][1][gene_dict[pmid][1].index(word) - 1] + " " + word)
        if stresslist:
            result_dict[pmid] = [gene_dict[pmid][0], stresslist]
        else:
            no_stress.append("PMID:" + pmid)

    return [result_dict, no_stress]


def search_co_oc(
        result_dict):  # Functie verantwoordelijk voor het maken van de uiteindelijke datastructuur die gebruikt word door de json filewriter
    new_result_dict = {}  # Input = result_dict,  output = completed_sunburst_dict
    dictionary = {}
    gene_stress_stringlist = []
    appearence_counter = []
    filtered_dict = {}
    final_dict = {}

    for PubMedID in result_dict:  # Deze for loop zorgt er voor dat een stress factor niet 2 keer in de stresslijst voor kan komen
        value = result_dict[PubMedID]
        genes = value[0]
        stresses = value[1]
        unique_stress_list = []
        for stress in stresses:
            if stress not in unique_stress_list:
                unique_stress_list.append(stress)
        new_result_dict[PubMedID] = [genes,
                                     unique_stress_list]  # dictionary key=pubmedID values= lijst van genen en lijst van unique stress omstandigheden

    for PubMedID in new_result_dict:  # Deze for loopt maakt in string in het formaat gen|-|stress deze string zal dienen als key voor een dictionary
        values = new_result_dict[PubMedID]  # een value uit new_result_dict bijv [[gen,gen1],[stres,stres1]]
        genes = values[0]
        stresses = values[1]
        for gene in genes:  # deze geneste for loop loopt over de gen lijst en stress lijst en maak een strings door gen en stress aan elkaar te koppelen middels een |-| teken
            for stress in stresses:
                gene_stress_stringlist_key = gene + "|-|" + stress
                gene_stress_stringlist.append(
                    gene_stress_stringlist_key)  # een lijst me de verbanden die te vinden zijn tussen gen en stress

    for key in gene_stress_stringlist:  # Deze for loopt maakt een dictionary in het formaat key=gen|-|stess value = hoevaak deze combinatie gevonden is
        if key not in appearence_counter:  # Als een key nog niet toegevoegt is aan dictionary
            appearence_counter.append(key)
            dictionary[key] = 1
        else:
            count = dictionary.get(key)  # Als de key al toegevoegd is aan dictionary word de value van deze opgehaald
            count += 1  # en word er 1 bij de value opgeteld
            dictionary[key] = count

    for key in dictionary:  # In deze for loopt word de uiteindelijke count van een combinatie gen|-|stress opgehaald uit de dictionary
        count = dictionary[key]

        for PubMedID in result_dict:  # In deze for loop word er een dictionary gemaakt in het formaat key=gen|-|stress value= hoevaak de combinatie voorkomt, in welke abstracts de combinatie voorkomt
            value = result_dict[PubMedID]
            genes = value[0]
            stresses = value[1]
            for gene in genes:
                for stress in stresses:
                    check = gene + "|-|" + stress  # Als de gen stress combinaite gelijk is aan de combinaite key in dictionary
                    if key == check:  # word het pubmed id(key)van result_dict toegevoegd aan een lijst dit word zo gedaan dat er alleen unique codes in de lijst aanwezig zijn
                        if key not in final_dict:
                            pubmed_id_list = [PubMedID]
                            final_dict[key] = count, pubmed_id_list
                        if key in final_dict:
                            acc_codes = final_dict.get(key)
                            acc_code = acc_codes[1]
                            if PubMedID not in acc_code:
                                acc_code.append(PubMedID)
                                final_dict[
                                    key] = count, acc_code  # Dictionary: key=gen|-|stress value= hoevaak de combinatie voorkomt, in welke abstracts de combinatie voorkomt
                            if PubMedID in acc_code:
                                final_dict[key] = count, acc_code

        for relation in final_dict:  # Filter om eventueel relaties weg te filteren als deze maar in 1 abstact voorkomen om zo significanter relaties over te houden
            appearances = final_dict[relation]
            if appearances[
                0] >= 1:  # Om het filter te uit te zetten gebruik >=, om het filter aan te zetten gebruik >, gewenste aantal voorkomen van een gen stress relatie
                filtered_dict[relation] = appearances

    stress_list = []
    sunburst_dict = {}

    for gene_stress in filtered_dict:  # Deze for loopt maakt een dic die geschikt is om een json file van te maken voor de sunburst. key=stress value=lijst van genen die een veband hebben met deze stress factor
        genes = gene_stress.split('|-|', 1)[0]
        stresses = gene_stress.split('|-|', 1)[1]
        stresses = stresses.lower()
        gene_list = []
        if stresses not in stress_list:  # Als stress nog niet in de stress_list voorkomt word deze een nieuwe key voor sunburst_dict
            for gene_stress_2 in filtered_dict:  # Deze loop voegd dictionarys samen afhankelijk van de stress
                genes_2 = gene_stress_2.split('|-|', 1)[0]
                stress_2 = gene_stress_2.split('|-|', 1)[1].lower()
                stress_list.append(stresses)
                if stresses == stress_2:  # Als de stresomstandigheid van 2 verbanden uit gen_stress hetzelfde zijn
                    if genes not in gene_list:  # Als gen 1 en 2 nog niet aanwezig zijn in genlist worden deze hieraan toegevoegd
                        gene_list.append(genes)
                    if genes_2 not in gene_list:
                        gene_list.append(genes_2)
            sunburst_dict[
                stresses] = gene_list  # Een dictionary. key= stressomstandigheid value=lijst van unique genen die een veband hebben met de stress omstandigheid

    completed_sunburst_dict = {}

    for stress_dict in sunburst_dict:  # deze for loopt maakt een completer dict die geschikt is om een json file van te maken voor de sunburst. Deze dictionary ziet er als als volgd uit {stress1{gen{[acc,acc,acc]}gen{[acc,acc,acc]}}stress2{gen{[acc,acc,acc]}}}
        gene_list = []  # deze dictionary word gemaakt door de accessiecodes van filterd_dict te koppelen aan een bijhorend gen in sunburst_dict.
        for gene_stress in filtered_dict:  # de koppeling vind plaats wanneer een het gen en stres uit gilterd_dict overeenkomen met het gen(value) en stress(key) uit sunburst_dict
            genes = gene_stress.split('|-|', 1)[0]
            stresses = gene_stress.split('|-|', 1)[1]
            if stresses == stress_dict:
                acc_list = filtered_dict[gene_stress][1]
                for word in sunburst_dict[stress_dict]:
                    gen_code_dict = {}
                    if word == genes:
                        gen_code_dict[word] = acc_list
                        gene_list.append(gen_code_dict)
        completed_sunburst_dict[stress_dict] = gene_list

    return completed_sunburst_dict  # Deze dictionary ziet er als als volgd uit {stress1{gen{[acc,acc,acc]}gen{[acc,acc,acc]}}stress2{gen{[acc,acc,acc]}}}


def json_file_writer(search, completed_sunburst_dict):
    file_handler = open('templates/flare.json', "w")
    file_handler.write("{\n")
    file_handler.write('"name": "' + str(search) + '", \n')
    file_handler.write('"children": [ \n')
    count = 0

    for stress in completed_sunburst_dict:
        count += 1
        file_handler.write("\t{\n")
        file_handler.write('\t"name": "' + stress + '",\n')
        file_handler.write('\t"children": [\n')
        index = 0

        for gene_dict in completed_sunburst_dict[stress]:
            for gene, acc_code in gene_dict.items():
                index += 1
                if index == len(completed_sunburst_dict[stress]):
                    file_handler.write('\t\t{"name": "' + gene + '", "size": ' + str(len(acc_code)) + '}\n')
                else:
                    file_handler.write('\t\t{"name": "' + gene + '", "size": ' + str(len(acc_code)) + '},\n')
        file_handler.write('\t\t]\n')
        if count == len(completed_sunburst_dict):
            file_handler.write('\t}\n')
        else:
            file_handler.write('\t},\n')

    file_handler.write('\t]\n')
    file_handler.write('}')
    file_handler.close()

    return "done"


def main(search):
    mail = "davidmaarn@gmail.com"
    query = search + " AND Anthocyanin"

    hits = get_number_of_hits(mail, query)
    records = get_records(mail, query, hits)
    info_dict = process_records(records)
    gene_dict = search_genes(info_dict[0])
    result_dict = search_stress(gene_dict[0])
    co_occurrence = search_co_oc(result_dict[0])
    done = json_file_writer(search, co_occurrence)

    return done


@app.route('/', methods=["GET"])
def mainpage():
    return render_template('index.html')


@app.route('/help', methods=["GET"])
def help_page():
    return render_template('help.html')


@app.route('/contact', methods=["GET"])
def contact_page():
    return render_template('contact.html')


@app.route('/visualisatie', methods=["GET"])
def visualisation_page():
    return render_template('visualisatie.html')


@app.route('/textmine', methods=['GET', 'POST'])
def textmine():
    if request.method == "POST":
        query = request.form['query']
        done = main(query)
        if done == "done":
            return render_template("textminedone.html", query=query)
    return render_template("textmine.html")


if __name__ == '__main__':
    app.run()