// C语言词法分析器
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <set>
#include <stack>
using namespace std;
/* 不要修改这个标准输入函数 */
void read_prog(string &prog)
{
    char c;
    while (scanf("%c", &c) != EOF)
    {
        prog += c;
    }
}
/* 你可以添加其他函数 */

enum class SymbolType
{
    TERMINAL,
    NON_TERMINAL,
    EPSILON,
    END
};

struct Symbol
{
    string name;
    SymbolType type;

    Symbol();
    Symbol(const string &n, SymbolType t);

    bool operator==(const Symbol &other) const;
    bool operator<(const Symbol &other) const;
};

struct Production
{
    Symbol left;          // A
    vector<Symbol> right; // α

    // A → α
    Production();
    Production(const Symbol &l, const vector<Symbol> &r);
};

class Grammar
{
public:
    // 构造
    Grammar();
    Grammar(const Symbol &start);

    void setStartSymbol(const Symbol &s);
    Symbol getStartSymbol() const;

    // 添加产生式 A → α
    void addProduction(const Symbol &left,
                       const vector<Symbol> &right);

    // 查询
    const vector<Production> &getProductions() const;
    const vector<Production> &getProductionsOf(const Symbol &nt) const;

    // 判断符号类型
    bool isTerminal(const Symbol &s) const;
    bool isNonTerminal(const Symbol &s) const;

    // 获取集合
    const set<Symbol> &getTerminals() const;
    const set<Symbol> &getNonTerminals() const;

private:
    Symbol startSymbol;
    vector<Production> productions;
    map<Symbol, vector<Production>> prodMap;

    set<Symbol> terminals;
    set<Symbol> nonTerminals;
};

class FirstFollowCalculator
{
public:
    FirstFollowCalculator(const Grammar &g);

    void computeFirst();
    void computeFollow();

    const set<Symbol> &getFirst(const Symbol &s) const;
    const set<Symbol> &getFollow(const Symbol &s) const;

    // FIRST(α)
    set<Symbol> firstOfSequence(const vector<Symbol> &seq) const;

private:
    const Grammar &grammar;

    map<Symbol, set<Symbol>> first;
    map<Symbol, set<Symbol>> follow;

    bool addToFirst(const Symbol &s, const Symbol &x);
    bool addToFollow(const Symbol &s, const Symbol &x);
};

class LL1Table
{
public:
    LL1Table();

    void build(const Grammar &g,
               const FirstFollowCalculator &ff);

    bool hasEntry(const Symbol &nonTerminal,
                  const Symbol &terminal) const;

    const Production &getEntry(const Symbol &nonTerminal,
                               const Symbol &terminal) const;

    bool hasConflict() const;

private:
    map<pair<Symbol, Symbol>, Production> table;
    bool conflict;
};

class LL1Parser
{
public:
    LL1Parser(const Grammar &g,
              const LL1Table &table);

    void parse(const vector<Symbol> &input);

private:
    const Grammar &grammar;
    const LL1Table &ll1Table;

    stack<Symbol> parseStack;

    void reportError(int line,
                     const string &msg);

    void synchronize(const Symbol &nonTerminal,
                     const vector<Symbol> &input,
                     size_t &pos);
};


void Analysis()
{
    string prog;
    read_prog(prog);
    /* 骚年们 请开始你们的表演 */
    /********* Begin *********/

    /********* End *********/
}