import os, re

with open(os.path.dirname(os.path.abspath(__file__)) + "/../output/push_test.txt", "r") as file:
  content = file.read()

assert all ( x not in content \
  for x in [ "warning", "Warning", "error", "Error", "Fail" ] ), \
  "Some tests seem to have failed!"

index = re.search("fail", content)
while index != None:
  assert index.start() > 8 and content[index.start()-9:index.start()] == ' 0 tests ', \
    "Some tests seem to have failed!"
  content = content[index.end():len(content)]
  index = re.search("fail", content)

print("All tests have passed!")
