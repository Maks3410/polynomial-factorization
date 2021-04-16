#include <vector>

using namespace std;

polynomial multiply_polynomials(polynomial a, polynomial b){
    polynomial ans{vector<int>(a.degree + b.degree + 1)};
    for (int i = 0; i <= a.degree; ++i) {
        for (int j = 0; j <= b.degree; ++j) {
            ans.monomials[i + j] += (a.monomials[i] * b.monomials[j]);
        }
    }
    return ans;
}

polynomial sum_polynomials(polynomial a, polynomial b) {
    polynomial answer{vector<int>(max(a.degree, b.degree) + 1)};
    int i;
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

success_and_result divide_poly_num(polynomial a, int n) {
    for (int i = 0; i <= a.degree; ++i) {
        if (a.monomials[i] % n != 0)
            return success_and_result{false};
        a.monomials[i] /= n;
    }
    return success_and_result{true, a};
}

success_and_result lagrange2(vector<vector<int>> points) {
    int n = points.size() - 1;
    polynomial answer{vector<int>(n)};
    for (int i = 0; i <= n; ++i) {
        polynomial li;
        int del = 1;
        for (int j = 0; j <= n; ++j) {
            if(i != j) {
                li = multiply_polynomials(li, polynomial{vector<int>{1, -points[j][0]}});
                del *= (points[i][0] - points[j][0]);
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
