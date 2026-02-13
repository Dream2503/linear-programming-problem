#pragma once

inline std::vector<std::vector<int>> lpp::detail::generate_combination(const int n, const int k) {
    std::vector<int> curr;
    std::vector<std::vector<int>> ans;

    const auto inner = [](auto&& self, int start, int n, int k, std::vector<int>& current, std::vector<std::vector<int>>& res) -> void {
        if (current.size() == k) {
            res.push_back(current);
            return;
        }

        for (int i = start; i < n; i++) {
            current.push_back(i);
            self(self, i + 1, n, k, current, res);
            current.pop_back();
        }
    };
    inner(inner, 0, n, k, curr, ans);
    return ans;
}
