from transformers import AutoTokenizer, AutoModelForTokenClassification

tokenizer = AutoTokenizer.from_pretrained("dslim/bert-base-NER")  # Example of a NER model
model = AutoModelForTokenClassification.from_pretrained("dslim/bert-base-NER")

text = "Mutations in the BRCA1 gene are associated with increased risk of breast cancer."
inputs = tokenizer(text, return_tensors="pt", truncation=True, padding=True)

outputs = model(**inputs)
predictions = outputs.logits.argmax(dim=-1)

# Example label mapping
label_map = {i: label for i, label in enumerate(tokenizer.get_vocab())}

# Decode predictions
predicted_labels = [label_map[p.item()] for p in predictions.squeeze()]

