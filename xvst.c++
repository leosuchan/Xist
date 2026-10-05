#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>
#include <queue>
#include <limits>
#include <numeric>
#include <algorithm>

namespace py = pybind11;
using namespace std;

struct EdgeInput {
    int u;
    int v;
    long long w;

    EdgeInput() : u(0), v(0), w(0) {}
    EdgeInput(int u_, int v_, long long w_)
        : u(u_), v(v_), w(w_) {}
};

struct XvstResult {
    double ncut_value;
    vector<vector<int>> partition;
};

struct DinicEdge {
    int to;
    int rev;
    long long cap;
};

class Dinic {
public:
    int N;
    vector<vector<DinicEdge>> G;
    vector<int> level;
    vector<int> ptr;
    vector<vector<long long>> init_caps;

    explicit Dinic(int n)
        : N(n), G(n), level(n), ptr(n) {}

    explicit Dinic(const vector<vector<DinicEdge>>& base)
        : N((int)base.size()), G(base), level(N), ptr(N) {
        init_caps.resize(N);
        for (int i = 0; i < N; ++i) {
            init_caps[i].resize(G[i].size());
            for (size_t j = 0; j < G[i].size(); ++j)
                init_caps[i][j] = G[i][j].cap;
        }
    }

    void restore_caps() {
        for (int i = 0; i < N; ++i)
            for (size_t j = 0; j < G[i].size(); ++j)
                G[i][j].cap = init_caps[i][j];
    }

    void add_edge(int u, int v, long long cap) {
        DinicEdge a{v, (int)G[v].size(), cap};
        DinicEdge b{u, (int)G[u].size(), 0};

        G[u].push_back(a);
        G[v].push_back(b);
    }

    bool bfs(int s, int t) {
        fill(level.begin(), level.end(), -1);

        queue<int> q;
        q.push(s);
        level[s] = 0;

        while (!q.empty()) {
            int v = q.front();
            q.pop();

            for (const auto& e : G[v]) {
                if (e.cap > 0 && level[e.to] == -1) {
                    level[e.to] = level[v] + 1;
                    q.push(e.to);
                }
            }
        }

        return level[t] != -1;
    }

    struct DfsFrame { int v; int edge; long long pushed; int parent_edge; };

    long long dfs(int v, int t, long long pushed) {
        if (pushed == 0) return 0;

        vector<DfsFrame> st;
        st.reserve(N);
        st.push_back({v, ptr[v], pushed, -1});

        while (!st.empty()) {
            DfsFrame &frame = st.back();
            int x = frame.v;

            if (x == t) {
                long long tr = frame.pushed;
                for (int i = 1; i < (int)st.size(); ++i) {
                    int u = st[i - 1].v;
                    int eidx = st[i].parent_edge;
                    DinicEdge &e = G[u][eidx];
                    e.cap -= tr;
                    G[e.to][e.rev].cap += tr;
                }
                return tr;
            }

            if (frame.edge >= (int)G[x].size()) {
                ptr[x] = frame.edge;
                st.pop_back();
                if (!st.empty()) {
                    st.back().edge += 1;
                }
                continue;
            }

            DinicEdge &e = G[x][frame.edge];
            if (level[e.to] != level[x] + 1 || e.cap <= 0) {
                frame.edge += 1;
                continue;
            }

            long long next_pushed = min(frame.pushed, e.cap);
            st.push_back({e.to, ptr[e.to], next_pushed, frame.edge});
        }

        return 0;
    }

    long long maxflow(int s, int t) {
        long long flow = 0;

        while (bfs(s, t)) {
            fill(ptr.begin(), ptr.end(), 0);
            while (long long pushed = dfs(s, t, numeric_limits<long long>::max())) {
                flow += pushed;
            }
        }

        return flow;
    }

    vector<bool> mincut_side(int s) {
        vector<bool> vis(N, false);
        queue<int> q;
        q.push(s);
        vis[s] = true;

        while (!q.empty()) {
            int u = q.front();
            q.pop();
            for (const auto& e : G[u]) {
                if (e.cap > 0 && !vis[e.to]) {
                    vis[e.to] = true;
                    q.push(e.to);
                }
            }
        }

        return vis;
    }
};

vector<long long> compute_strengths(
    int n,
    const vector<EdgeInput>& edges
) {
    vector<long long> deg(n, 0);
    for (const auto& e : edges) {
        deg[e.u] += e.w;
        deg[e.v] += e.w;
    }
    return deg;
}



XvstResult xvst2(
    int n,
    const vector<EdgeInput>& edges
) {
    double best_ncut = numeric_limits<double>::infinity();
    vector<vector<int>> best_partition;

    if (n <= 1) {
        return {best_ncut, best_partition};
    }

    vector<long long> deg = compute_strengths(n, edges);
    long long total_volume = accumulate(deg.begin(), deg.end(), 0LL);

    vector<int> incident(n, 0);
    for (const auto &e : edges) {
        incident[e.u]++;
        incident[e.v]++;
    }

    vector<vector<DinicEdge>> baseG(n);
    for (int i = 0; i < n; ++i) {
        baseG[i].reserve(incident[i] * 2 + 2);
    }

    auto add_to_base = [&](int u, int v, long long cap) {
        DinicEdge a{v, (int)baseG[v].size(), cap};
        DinicEdge b{u, (int)baseG[u].size(), 0};
        baseG[u].push_back(a);
        baseG[v].push_back(b);
    };

    for (const auto &e : edges) {
        add_to_base(e.u, e.v, e.w);
        add_to_base(e.v, e.u, e.w);
    }

    Dinic dinic(baseG);
    for (int i = 0; i < n - 1; ++i) {
        for (int j = i + 1; j < n; ++j) {
            dinic.restore_caps();
            long long flow = dinic.maxflow(i, j);
            vector<bool> visited = dinic.mincut_side(i);

            vector<int> partA, partB;
            partA.reserve(n);
            partB.reserve(n);
            for (int v = 0; v < n; ++v) {
                if (visited[v]) {
                    partA.push_back(v);
                } else {
                    partB.push_back(v);
                }
            }

            const vector<int>& part_i = visited[i] ? partA : partB;

            long long vol_part_i = 0;
            for (int v : part_i) {
                vol_part_i += deg[v];
            }
            long long vol_part_i_complement = total_volume - vol_part_i;

            if (vol_part_i > 0 && vol_part_i_complement > 0) {
                double ncut = static_cast<double>(flow) /
                    (static_cast<double>(vol_part_i) * static_cast<double>(vol_part_i_complement));
                if (ncut < best_ncut) {
                    best_ncut = ncut;
                    best_partition = {partA, partB};
                }
            }
        }
    }

    return {best_ncut, best_partition};
}

py::list xvst_from_edges2(
    int n,
    const vector<EdgeInput>& edges
) {
    XvstResult result = xvst2(n, edges);
    py::list output;
    output.append(result.ncut_value);
    output.append(result.partition);
    return output;
}

PYBIND11_MODULE(xvst_dinic_faster, m) {
    py::class_<EdgeInput>(m, "EdgeInput", py::module_local())
        .def(py::init<>())
        .def(py::init<int, int, long long>())
        .def_readwrite("u", &EdgeInput::u)
        .def_readwrite("v", &EdgeInput::v)
        .def_readwrite("w", &EdgeInput::w);

    m.def(
        "xvst",
        py::overload_cast<
            int,
            const vector<EdgeInput>&
        >(&xvst_from_edges2)
    );
}