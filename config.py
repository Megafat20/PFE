OPEN_API="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"
key="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"

miralAI="eVNVoFHuxr0HU2p9RtRnMosg6q1uibQN"
m2="DUjNVQ9Lm92jk76pEB6UsFjJh9qYPEIx"

llama_api="09b970c5-aaff-4f92-be81-d50506f3153c"
llama_cloud_api="llx-6Bv33PpqkpYKpAsKuH2ecVW2dAZZ1ukMl9EGa7YfQRJhHabU"


"""  
# Configure Mistral with rate limiting
Settings.llm = MistralAI(
    temperature=0.8,
    model="mistral-small-latest",  # Start with smaller model for testing
    max_retries=3,
    timeout=30,
    api_key="DUjNVQ9Lm92jk76pEB6UsFjJh9qYPEIx"
)

# Configure local embeddings
Settings.embed_model = "local:BAAI/bge-small-en-v1.5"
Settings.chunk_size = 512  # Reduce processing load


from openai import OpenAI

client = OpenAI(
  base_url="https://openrouter.ai/api/v1",
  api_key="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",
)

completion = client.chat.completions.create(
  extra_headers={
    "HTTP-Referer": "<YOUR_SITE_URL>", # Optional. Site URL for rankings on openrouter.ai.
    "X-Title": "<YOUR_SITE_NAME>", # Optional. Site title for rankings on openrouter.ai.
  },
  extra_body={},
  model="deepseek/deepseek-r1:free",
  messages=[
    {
      "role": "user",
      "content": "What is the meaning of life?"
    }
  ]
)
print(completion.choices[0].message.content)

"""