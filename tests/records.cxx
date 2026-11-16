// tests/records.hxx: parse and compare golden records in-process, one case per rule.
#include <tests/records.hxx>

#include <cstdio>
#include <string>
#include <vector>

int failures = 0;
void check(const bool ok, const char* what)
{
  std::printf("%s  %s\n", ok ? "ok  " : "FAIL", what);
  failures += !ok;
}

// The comparison of the printed lines with the expected text at the default tolerances, its
// complaints printed for the log.
tests::Comparison compare(const std::vector<std::string>& printed, const std::string& expected)
{
  const tests::Comparison result = tests::compare(printed, expected, 1e-5, 0.1);
  for (const std::string& complaint : result.complaints)
    std::printf("        complaint: %s\n", complaint.c_str());
  return result;
}
bool mentions(const tests::Comparison& result, const char* text)
{
  for (const std::string& complaint : result.complaints)
    if (complaint.find(text) != std::string::npos)
      return true;
  return false;
}

int main()
{
  const std::vector<std::string> printed = {
    "- {dim: 2, degree: 1, n: 2, dofs: 24, cg: converged_rtol, iterations: 5, error: "
    "1.107770e-02}\n",
    "- {dim: 2, degree: 1, n: 4, dofs: 80, cg: converged_rtol, iterations: 21, error: "
    "3.553054e-03}\n"};

  tests::Record record;
  check(tests::parse_record(printed[0], record) && record.size() == 7 &&
          record[4] == std::make_pair(std::string("cg"), std::string("converged_rtol")),
        "parse a record line");
  check(tests::to_string(record) + "\n" == printed[0], "to_string gives the printed line back");
  check(!tests::parse_record("  direct vs cg: 1.2e-15", record), "other lines are no records");

  {
    const auto result = compare(printed, printed[0] + printed[1]);
    check(result.complaints.empty() && result.checked == 2 && result.unchecked == 0,
          "the printed records themselves pass");
  }
  {
    // a seed: fewer keys and records, digits of another implementation
    const auto result =
      compare(printed, "- {dim: 2, degree: 1, n: 4, error: 0.0035530543400028873}\r\n");
    check(result.complaints.empty() && result.checked == 1 && result.unchecked == 1,
          "a seed with fewer keys and records, CRLF");
  }
  {
    const auto result = compare(printed, "- {dim: 2, degree: 1, n: 4, error: 3.6e-03}");
    check(result.complaints.size() == 1 && mentions(result, "n: 4") &&
            mentions(result, "error is 3.553054e-03, expected 3.6e-03"),
          "a float beyond rtol names the record and the key");
  }
  {
    const auto result = compare(printed, "- {dim: 2, degree: 1, n: 2, cg: converged_atol}");
    check(result.complaints.size() == 1 &&
            mentions(result, "cg is converged_rtol, expected converged_atol"),
          "a changed stop reason: no record matches, the nearest names the key");
  }
  {
    const auto within = compare(printed, "- {n: 4, iterations: 23}");
    const auto beyond = compare(printed, "- {n: 4, iterations: 26}");
    check(within.complaints.empty() && beyond.complaints.size() == 1 &&
            mentions(beyond, "iterations is 21, expected 26"),
          "iterations: within 0.1 |e| + 2 passes, beyond fails");
  }
  {
    const auto result = compare(printed, "- {n: 4, residual: 1e-13}");
    check(result.complaints.size() == 1 && mentions(result, "residual missing"),
          "a compared key missing from the printed record");
  }
  {
    const auto result = compare(printed, "- {n: 8, error: 1e-3}");
    check(result.complaints.size() == 1 && mentions(result, "n is 2, expected 8"),
          "an expected record that was not printed");
  }
  {
    const auto result = compare(printed, "- {dim: 2, error: 1e-3}");
    check(result.complaints.size() == 1 && mentions(result, "2 printed records match"),
          "an ambiguous expected record");
  }
  {
    const auto result = compare(printed, "# nothing\n");
    check(result.complaints.size() == 1 && mentions(result, "holds no records"),
          "an expected file without records");
  }
  {
    const auto result = compare({}, "- {n: 2, error: 1e-2}");
    check(result.complaints.size() == 1 && mentions(result, "no record printed"),
          "no record printed");
  }
  {
    const auto result = compare({"- {n: 2, error: nan}"}, "- {n: 2, error: 1e-2}");
    check(result.complaints.size() == 1, "nan fails");
  }
  return failures;
}
