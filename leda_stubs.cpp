#include <LEDA/numbers/integer.h>
#include <LEDA/numbers/rational.h>
#include <LEDA/numbers/real.h>
#include <LEDA/numbers/vector.h>
#include <LEDA/core/string.h>

int leda_set_fpu_defaults() { return 0; }

namespace leda {
    int integer::cmp(const integer&, const integer&) { return 0; }
    int rational::cmp(const rational&, const rational&) { return 0; }
    int compare(const real&, const real&) { return 0; }
    vector::vector(int d) {}
    bool is_file(string) { return false; }
    bool delete_file(string) { return false; }
    string tmp_file_name() { return string(""); }
}
