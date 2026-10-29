import os, re

with open(os.path.dirname(os.path.abspath(__file__)) + "/../output/reproducibles_test.txt", "r") \
  as file:
  content = file.read()

## Messages that indicate a failed run. Note that a message like "CG solver failed ... Trying
#  GMRES" is not a failure as long as the fallback solver converges (otherwise "also failed" or
#  "ERROR" is printed).
FAILURE_PATTERNS = [ "Traceback", "Assert failed", "ERROR", "also failed", "did not converge",
                     "warning:", "Warning:", "error:" ]

failed = False
for script in re.split(r"\n\nTry ", content):
  name = script.split(" ...", 1)[0].strip()
  problems = [ pattern for pattern in FAILURE_PATTERNS if pattern in script ]
  # "Error" must only be used for printed errors of a discretization, i.e., followed by ': '.
  if any(script[match.end():match.end() + 2] != ': ' for match in re.finditer("Error", script)):
    problems.append("Error")
  if problems:
    failed = True
    print("FAILED:", os.path.basename(name), "(", ", ".join(problems), ")")

assert not failed, "Some tests seem to have failed!"
print("All tests have passed!")
