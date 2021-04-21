import os, re

with open(os.path.dirname(os.path.abspath(__file__)) + "/../output/reproducibles_test.txt", "r") \
  as file:
  content = file.read()

assert ("warning" or "Warning" or "error" or "fail" or "Fail" or "assert" or "Assert") \
  not in content, "Some tests seem to have failed!"

index = re.search("Error", content)
while index != None:
  assert content[index.end():index.end() + 2] == ': ', "Some tests seem to have failed!"
  content = content[index.end():len(content)]
  index = re.search("Error", content)

print("All tests have passed!")
