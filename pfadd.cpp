#include <vector>

using namespace std;

polynomial multiply_polynomials(polynomial a, polynomial b){
    polynomial ans{vector<long long>(a.degree + b.degree + 1)};
    for (int i = 0; i <= a.degree; ++i) {
        for (int j = 0; j <= b.degree; ++j) {
            ans.monomials[i + j] += (a.monomials[i] * b.monomials[j]);
        }
    }
    return ans;
}

polynomial sum_polynomials(polynomial a, polynomial b) {
    polynomial answer{vector<long long>(max(a.degree, b.degree) + 1)};
    long long i;
    for (i = 0; i <= min(a.degree, b.degree); ++i) {
        answer.monomials[i] = a.monomials[i] + b.monomials[i];
    }
    if (a.degree > b.degree) {
        while (i <= answer.degree) {
            answer.monomials[i] = a.monomials[i];
        }
    }
    else {
        while (i <= answer.degree) {
            answer.monomials[i] = b.monomials[i];
        }
    }
    return answer;
}

success_and_result divide_poly_num(polynomial a, long long n) {
    for (long long i = 0; i <= a.degree; ++i) {
        if (a.monomials[i] % n != 0)
            return success_and_result{false};
        a.monomials[i] /= n;
    }
    return success_and_result{true, a};
}

success_and_result lagrange2(vector<long long> x, vector<long long> y) {
    long long n = x.size() - 1;
    // polynomial li{vector<long long>(n)};
    polynomial answer;
    for (int i = 0; i <= n; ++i) {
        polynomial li{vector<long long>(n)};
        long long del = 1;
        for (int j = 0; j <= n; ++j) {
            if(i != j) {
                li = multiply_polynomials(li, polynomial{vector<long long>{1, -x[j]}});
                del *= (x[i] - x[j]);
            }
        }
        auto division = divide_poly_num(li, del);
        if (division.success) {
            li = division.result;
            answer = sum_polynomials(answer, li);
        }
        else
            return success_and_result{false};
    }
    return {true, answer};
}

template <typename T>
void print_vec(const vector<T>& a){
    for (auto i : a) {
        cout << i << ' ';
    }
    cout << endl;
}
