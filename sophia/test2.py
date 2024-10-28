import spacy
from collections import Counter

# Load spaCy model
nlp = spacy.load("en_core_web_sm")


# Your text here
text = """
TP53 mutations are commonly associated with various cancers including breast and lung cancer.
The BRCA1 gene plays a critical role in DNA repair, and mutations in this gene can lead to an increased risk of breast and ovarian cancer.
Mutations in the EGFR gene are often observed in non-small cell lung cancer and can impact treatment responses.
BRCA2 mutations have been linked to a higher risk of prostate and pancreatic cancers.
In lung cancer, EGFR mutations can be targeted by specific inhibitors for treatment.
"""

# Process the text
doc = nlp(text)

# Extract entities and their labels
entities = [(ent.text, ent.label_) for ent in doc.ents]
print(entities)

# Count entity frequencies
entity_counts = Counter(entities)

# Get the most common entities
top_entities = entity_counts.most_common(10)

# Separate entities and counts for plotting
entity_labels = [entity[0] for entity in top_entities]
