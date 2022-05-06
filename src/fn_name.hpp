#ifndef FN_NAME_HPP
#define FN_NAME_HPP

#include <algorithm>
#include <iostream>
#include <iterator>
#include <stdexcept>

template<typename InputIterator, typename T>
InputIterator
find_closing(InputIterator first, InputIterator last, T close)
{
  if (first == last) {
    return last;
  }

  auto open = *first;
  unsigned counter = 1;
  while (++first != last) {
    if (*first == close && --counter == 0) {
      return first;
    }
    if (*first == open) {
      ++counter;
    }
  }

  return last;
}

template<std::size_t N, std::size_t N2>
std::string
get_fn_name(char const (&str)[N], char const (&name)[N2])
{
  using std::begin, std::end, std::invalid_argument, std::reverse_iterator,
    std::search;

  // Argument to isalnum must be unsigned:
  auto cond = [](unsigned char c) { return !isalnum(c) && c != '_'; };

  auto iter = str;
  for (;; ++iter) {
    iter = search(iter, end(str), begin(name), end(name) - 1);

    if (iter == end(str)) {
      throw invalid_argument("");
    }

    if ((iter == begin(str) || cond(iter[-1])) &&
        (iter == end(str) - N2 ||
         (cond(iter[N2 - 1]) && iter[N2 - 1] != ':'))) {
      break;
    }
  }

  auto origin_iter = iter;
  while (iter != begin(str)) {
    --iter;
    for (const auto* p : { "()", "{}" }) {
      if (*iter == p[1]) {
        iter = find_closing(reverse_iterator<char const*>(iter + 1),
                            reverse_iterator<char const*>(begin(str)),
                            p[0])
                 .base() -
               2;
      }
    }

    if (cond(*iter) && *iter != ':') {
      return std::string(iter + 1, origin_iter + N2 - 1);
    }
  }

  return std::string(iter, origin_iter + N2 - 1);
}

#define FN_NAME get_fn_name(__PRETTY_FUNCTION__, __func__)

#endif