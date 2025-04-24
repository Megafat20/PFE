from llama_index.core import VectorStoreIndex, SimpleDirectoryReader
from llama_index.core.query_engine import CustomQueryEngine
from llama_index.core.prompts import PromptTemplate

# Load and index documents
docs = SimpleDirectoryReader(r"C:\Users\LENOVO\Documents\My Web Sites\Chatbot\data\arxiv\2310.07282v2.pdf").load_data()
index = VectorStoreIndex.from_documents(docs)

# Define an opinion-based prompt
opinion_template = PromptTemplate(
    "You are an expert in this domain. Give your opinion on the following question based on the document context:\n\n{context_str}\n\nQuestion: {query_str}\n\nOpinion:"
)

# Create a custom query engine
query_engine = index.as_query_engine(text_qa_template=opinion_template)

# Ask an opinion-based question
response = query_engine.query("What is your opinion on the approach used in this research paper?")
print(response)
