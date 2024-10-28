import spacy
from pathlib import Path
from collections import Counter
import matplotlib.pyplot as plt

from transformers import pipeline

# Load BioBERT model from HuggingFace for NER
nlp = pipeline("ner", model="dmis-lab/biobert-base-cased-v1.1")

# Load spaCy model
# nlp = spacy.load("en_core_sci_sm")

# Load the text from the files into the variable "text"
file_paths = [
    '/Users/sq566/Desktop/sophia/article_1.txt',
    '/Users/sq566/Desktop/sophia/article_2.txt',
    '/Users/sq566/Desktop/sophia/article_3.txt',
    '/Users/sq566/Desktop/sophia/article_4.txt',
    '/Users/sq566/Desktop/sophia/article_5.txt'
]

# Read and combine all the text from the files into one variable
text = "\n".join(Path(file).read_text() for file in file_paths)

# Process the text
doc = nlp(text)

# Extract entities and their labels
entities = [ent.text for ent in doc.ents]

# Count entity frequencies
entity_counts = Counter(entities)

# Get the most common entities
top_entities = entity_counts.most_common(10)

# Separate entities and counts for plotting
entity_labels = [entity[0] for entity in top_entities]
entity_frequencies = [count for entity, count in top_entities]

# Create the plot for entity labels and their frequencies
plt.figure(figsize=(10, 6))
plt.bar(entity_labels, entity_frequencies, color='skyblue')
plt.ylabel('Count')
plt.xlabel('Entities')
plt.title('Common biological entities')
plt.show()

