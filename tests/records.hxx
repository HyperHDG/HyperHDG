#pragma once  // Ensure that file is included only once in a single compilation.

/*!*************************************************************************************************
 * \file    tests/records.hxx
 * \brief   Golden records: the record lines a test prints, checked against an expected file.
 *
 * A record is one line `- {key: value, ...}`, printed by tests::record with a printf format. At the
 * end, tests::check compares the records with RECORDS_FILE (tests/expected/<ctest name>.yaml, set
 * by tests/CMakeLists.txt for the tests registered with RECORDS):
 * - An expected record matches the one printed record that agrees with it in every integer and
 *   string key it has (except iterations); floats are compared to rtol, iterations to
 *   iterations_rtol |e| + 2. Strings such as the stop reason of CG are exact.
 * - An expected record may hold fewer keys than the printed one, and printed records without an
 *   expected one are only counted, so a file can pin a part of the output (the seeds from the
 *   Cython-era scripts hold {dim, degree, n, error}).
 * - With HYPERHDG_RECORD=1 in the environment, check writes the printed records to RECORDS_FILE
 *   instead: `HYPERHDG_RECORD=1 ctest -R timoshenko` re-records a set.
 **************************************************************************************************/

#include <algorithm>
#include <cmath>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace tests
{
/*!*************************************************************************************************
 * \brief   A parsed record: its keys and values as written, in order.
 **************************************************************************************************/
using Record = std::vector<std::pair<std::string, std::string> >;

/*!*************************************************************************************************
 * \brief   The record of a line `- {key: value, ...}`, or false if the line is none.
 **************************************************************************************************/
inline bool parse_record(std::string line, Record& record)
{
  while (!line.empty() && (line.back() == '\n' || line.back() == '\r' || line.back() == ' '))
    line.pop_back();
  if (line.compare(0, 3, "- {") != 0 || line.back() != '}')
    return false;
  record.clear();
  std::stringstream fields(line.substr(3, line.size() - 4));
  for (std::string field; std::getline(fields, field, ',');)
  {
    const std::size_t begin = field.find_first_not_of(' '), colon = field.find(": ");
    if (begin == std::string::npos || colon == std::string::npos)
      return false;
    record.emplace_back(field.substr(begin, colon - begin), field.substr(colon + 2));
  }
  return !record.empty();
}

/*!*************************************************************************************************
 * \brief   The line of a record, as the test printed it (for messages and the expected file).
 **************************************************************************************************/
inline std::string to_string(const Record& record)
{
  std::string line = "- {";
  for (const auto& [key, value] : record)
    line += (line.size() > 3 ? ", " : "") + key + ": " + value;
  return line + "}";
}

namespace detail
{
enum class Kind
{
  integer,
  floating,
  string
};
inline Kind kind(const std::string& value)
{
  char* end = nullptr;
  std::strtol(value.c_str(), &end, 10);
  if (!value.empty() && *end == '\0')
    return Kind::integer;
  std::strtod(value.c_str(), &end);
  if (!value.empty() && *end == '\0')
    return Kind::floating;
  return Kind::string;
}
// The keys that identify a record: integers and strings, but not iterations.
inline bool identifies(const std::string& key, const std::string& value)
{
  return key != "iterations" && kind(value) != Kind::floating;
}
inline const std::string* find(const Record& record, const std::string& key)
{
  for (const auto& [k, v] : record)
    if (k == key)
      return &v;
  return nullptr;
}
inline bool same(const std::string& a, const std::string& b)
{
  return kind(a) == Kind::integer && kind(b) == Kind::integer
           ? std::strtoll(a.c_str(), nullptr, 10) == std::strtoll(b.c_str(), nullptr, 10)
           : a == b;
}
// The identity keys of expected in which printed differs (missing or another value).
inline std::vector<std::string> differences(const Record& expected, const Record& printed)
{
  std::vector<std::string> out;
  for (const auto& [key, value] : expected)
  {
    if (!identifies(key, value))
      continue;
    const std::string* other = find(printed, key);
    if (!other)
      out.push_back(key + " missing");
    else if (!same(*other, value))
      out.push_back(key + " is " + *other + ", expected " + value);
  }
  return out;
}
inline std::vector<std::string>& printed_lines()
{
  static std::vector<std::string> lines;
  return lines;
}
}  // end of namespace detail

/*!*************************************************************************************************
 * \brief   What compare found: the complaints, and how many printed records were checked or not.
 **************************************************************************************************/
struct Comparison
{
  std::vector<std::string> complaints;
  unsigned int checked = 0, unchecked = 0;
};

/*!*************************************************************************************************
 * \brief   Compare printed record lines with the text of an expected file (rules: see the top).
 **************************************************************************************************/
inline Comparison compare(const std::vector<std::string>& printed_lines,
                          const std::string& expected_text,
                          const double rtol,
                          const double iterations_rtol)
{
  Comparison result;
  std::vector<Record> printed, expected;
  Record record;
  for (const std::string& line : printed_lines)
    if (parse_record(line, record))
      printed.push_back(record);
  std::stringstream text(expected_text);
  for (std::string line; std::getline(text, line);)
    if (parse_record(line, record))
      expected.push_back(record);
  if (expected.empty())
    result.complaints.push_back("the expected file holds no records");

  std::vector<bool> matched(printed.size(), false);
  for (const Record& e : expected)
  {
    std::vector<std::size_t> candidates;
    for (std::size_t k = 0; k < printed.size(); ++k)
      if (detail::differences(e, printed[k]).empty())
        candidates.push_back(k);
    if (candidates.size() != 1)
    {
      std::string complaint = to_string(e) + ": ";
      if (candidates.size() > 1)
        complaint += std::to_string(candidates.size()) + " printed records match";
      else if (printed.empty())
        complaint += "no record printed";
      else
      {
        // the nearest printed record: the fewest differing identity keys
        std::vector<std::string> nearest;
        for (const Record& p : printed)
        {
          const std::vector<std::string> d = detail::differences(e, p);
          if (nearest.empty() || d.size() < nearest.size())
            nearest = d;
        }
        complaint += "no printed record matches; the nearest:";
        for (const std::string& d : nearest)
          complaint += (&d == &nearest.front() ? " " : "; ") + d;
      }
      result.complaints.push_back(complaint);
      continue;
    }

    const Record& p = printed[candidates.front()];
    matched[candidates.front()] = true;
    ++result.checked;
    for (const auto& [key, value] : e)
    {
      if (detail::identifies(key, value))
        continue;
      const std::string* other = detail::find(p, key);
      if (!other)
      {
        result.complaints.push_back(to_string(e) + ": " + key + " missing");
        continue;
      }
      const double a = std::strtod(other->c_str(), nullptr),
                   x = std::strtod(value.c_str(), nullptr);
      const double tolerance =
        key == "iterations" ? iterations_rtol * std::abs(x) + 2. : rtol * std::abs(x);
      if (!(std::abs(a - x) <= tolerance))  // also catches nan
        result.complaints.push_back(to_string(e) + ": " + key + " is " + *other + ", expected " +
                                    value);
    }
  }
  result.unchecked = std::count(matched.begin(), matched.end(), false);
  return result;
}

/*!*************************************************************************************************
 * \brief   Print a record line (printf format) and keep it for check.
 **************************************************************************************************/
[[gnu::format(printf, 1, 2)]] inline void record(const char* format, ...)
{
  std::va_list args, copy;
  va_start(args, format);
  va_copy(copy, args);
  const int length = std::vsnprintf(nullptr, 0, format, copy);
  va_end(copy);
  std::vector<char> line(std::max(length, 0) + 1, '\0');
  std::vsnprintf(line.data(), line.size(), format, args);
  va_end(args);
  std::fputs(line.data(), stdout);
  detail::printed_lines().emplace_back(line.data());
}

/*!*************************************************************************************************
 * \brief   Check the kept records against RECORDS_FILE, or write it with HYPERHDG_RECORD=1; returns
 *          failures plus the number of complaints. Without RECORDS_FILE (a compile outside CMake)
 *          nothing is checked.
 **************************************************************************************************/
inline int check(const int failures, const double rtol = 1e-5, const double iterations_rtol = 0.1)
{
#ifndef RECORDS_FILE
  (void)rtol;
  (void)iterations_rtol;
  std::printf("records: not checked (no RECORDS_FILE)\n");
  return failures;
#else
  const std::vector<std::string>& lines = detail::printed_lines();
  const char* record_mode = std::getenv("HYPERHDG_RECORD");
  if (record_mode && std::string(record_mode) == "1")
  {
    if (failures > 0 || lines.empty())
    {
      std::printf("records: not written to %s (%d failures, %zu records)\n", RECORDS_FILE, failures,
                  lines.size());
      return failures + 1;
    }
    std::ofstream file(RECORDS_FILE);
    Record record;
    for (const std::string& line : lines)
      if (parse_record(line, record))
        file << to_string(record) << '\n';
    if (!file)
    {
      std::printf("records: cannot write %s\n", RECORDS_FILE);
      return 1;
    }
    std::printf("records: %zu written to %s\n", lines.size(), RECORDS_FILE);
    return 0;
  }

  std::ifstream file(RECORDS_FILE);
  if (!file)
  {
    std::printf("records: cannot read %s (record it with HYPERHDG_RECORD=1)\n", RECORDS_FILE);
    return failures + 1;
  }
  std::stringstream text;
  text << file.rdbuf();
  const Comparison result = compare(lines, text.str(), rtol, iterations_rtol);
  for (const std::string& complaint : result.complaints)
    std::printf("  FAIL: %s\n", complaint.c_str());
  std::printf("records: %u checked, %u unchecked, %zu complaints\n", result.checked,
              result.unchecked, result.complaints.size());
  // an exit status counts modulo 256: keep the sum small
  return std::min(failures + (int)result.complaints.size(), 100);
#endif
}
}  // end of namespace tests
