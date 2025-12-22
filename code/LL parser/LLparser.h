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

// ========= enum =========

enum class SymbolType
{
    TERMINAL,
    NON_TERMINAL,
    EPSILON,
    END
};

// ========= Symbol =========

struct Symbol
{
    string name;
    SymbolType type;

    Symbol();
    Symbol(const string &n, SymbolType t);

    bool operator==(const Symbol &other) const;
    bool operator!=(const Symbol &other) const;
    bool operator<(const Symbol &other) const;
};

// ========= 规范ε和$符号 =========

const Symbol EPSILON("ε", SymbolType::EPSILON);
const Symbol END_MARK("$", SymbolType::END);

// ========= Production =========

struct Production
{
    Symbol left;          // A
    vector<Symbol> right; // α

    // A → α
    Production();
    Production(const Symbol &l, const vector<Symbol> &r);
};

// ========= Grammar =========

class Grammar
{
public:
    // 构造
    Grammar();
    Grammar(const Symbol &start);

    void setStartSymbol(const Symbol &s);
    Symbol getStartSymbol() const;

    // 添加产生式 A → α
    void addProduction(const Symbol &left, const vector<Symbol> &right);

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

// ========= FirstFollowCalculator =========

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

// ========= LL1Table =========

class LL1Table
{
public:
    LL1Table();

    void build(const Grammar &g, const FirstFollowCalculator &ff);

    bool hasEntry(const Symbol &nonTerminal, const Symbol &terminal) const;
    const Production &getEntry(const Symbol &nonTerminal, const Symbol &terminal) const;

    bool hasConflict() const;

private:
    map<pair<Symbol, Symbol>, Production> table;
    bool conflict;
};

// ========= LL1Parser =========

class LL1Parser
{
public:
    LL1Parser(const Grammar &g, const LL1Table &table);

    void parse(const vector<Symbol> &input);

private:
    const Grammar &grammar;
    const LL1Table &ll1Table;

    stack<Symbol> parseStack;

    void reportError(int line, const string &msg);

    void synchronize(const Symbol &nonTerminal, const vector<Symbol> &input,
                     size_t &pos);
};

// ========= 实现 =========

Symbol::Symbol() : name(""), type(SymbolType::TERMINAL) {}

Symbol::Symbol(const string &n, SymbolType t) : name(n), type(t) {}

bool Symbol::operator==(const Symbol &other) const
{
    return name == other.name && type == other.type;
}

bool Symbol::operator!=(const Symbol &other) const
{
    return name != other.name && type != other.type;
}

bool Symbol::operator<(const Symbol &other) const
{
    if (type != other.type)
        return type < other.type;
    return name < other.name;
}

// ========= 实现 =========

Production::Production() {}

Production::Production(const Symbol &l, const vector<Symbol> &r) : left(l), right(r) {}

// ========= 实现 =========

Grammar::Grammar() {}

Grammar::Grammar(const Symbol &start) : startSymbol(start)
{
    nonTerminals.insert(start);
}

void Grammar::setStartSymbol(const Symbol &s)
{
    startSymbol = s;
    nonTerminals.insert(s);
}

Symbol Grammar::getStartSymbol() const
{
    return startSymbol;
}

void Grammar::addProduction(const Symbol &left, const vector<Symbol> &right)
{
    // 左部一定是非终结符
    nonTerminals.insert(left);

    Production p(left, right);
    productions.push_back(p);
    prodMap[left].push_back(p);

    // 处理右部符号
    for (const Symbol &s : right)
    {
        if (s.type == SymbolType::NON_TERMINAL)
        {
            nonTerminals.insert(s);
        }
        else if (s.type == SymbolType::TERMINAL)
        {
            terminals.insert(s);
        }
        // EPSILON 不放入 terminals
    }
}

const vector<Production> &Grammar::getProductions() const
{
    return productions;
}

const vector<Production> &Grammar::getProductionsOf(const Symbol &nt) const
{
    static vector<Production> empty;
    auto it = prodMap.find(nt);
    if (it == prodMap.end())
        return empty;
    return it->second;
}

bool Grammar::isTerminal(const Symbol &s) const
{
    return terminals.count(s) > 0;
}

bool Grammar::isNonTerminal(const Symbol &s) const
{
    return nonTerminals.count(s) > 0;
}

const set<Symbol> &Grammar::getTerminals() const
{
    return terminals;
}

const set<Symbol> &Grammar::getNonTerminals() const
{
    return nonTerminals;
}

// ========= 实现 =========

FirstFollowCalculator::FirstFollowCalculator(const Grammar &g) : grammar(g)
{
    // 初始化 FIRST / FOLLOW 表
    for (const Symbol &nt : grammar.getNonTerminals())
    {
        first[nt] = set<Symbol>();
        follow[nt] = set<Symbol>();
    }

    // 终结符的 FIRST 是它自己
    for (const Symbol &t : grammar.getTerminals())
    {
        first[t].insert(t);
    }
}

const set<Symbol> &FirstFollowCalculator::getFirst(const Symbol &s) const
{
    static set<Symbol> empty;
    auto it = first.find(s);
    if (it == first.end())
        return empty;
    return it->second;
}

const set<Symbol> &FirstFollowCalculator::getFollow(const Symbol &s) const
{
    static set<Symbol> empty;
    auto it = follow.find(s);
    if (it == follow.end())
        return empty;
    return it->second;
}

set<Symbol> FirstFollowCalculator::firstOfSequence(const vector<Symbol> &seq) const
{
    set<Symbol> result;

    // 空串 ⇒ ε
    if (seq.empty())
    {
        result.insert(EPSILON);
        return result;
    }

    bool allNullable = true;

    for (const Symbol &s : seq)
    {
        const set<Symbol> &fs = getFirst(s);

        for (const Symbol &x : fs)
        {
            if (x != EPSILON)
                result.insert(x);
        }

        if (fs.count(EPSILON) == 0)
        {
            allNullable = false;
            break;
        }
    }

    if (allNullable)
    {
        result.insert(EPSILON);
    }

    return result;
}

void FirstFollowCalculator::computeFirst()
{
    bool changed = true;

    // EPSILON 的 FIRST 是它自己
    first[EPSILON].insert(EPSILON);

    while (changed)
    {
        changed = false;

        for (const Production &p : grammar.getProductions())
        {
            const Symbol &A = p.left;
            const vector<Symbol> &alpha = p.right;

            bool allNullable = true;

            for (const Symbol &X : alpha)
            {
                for (const Symbol &s : getFirst(X))
                {
                    if (s != EPSILON)
                    {
                        if (first[A].insert(s).second)
                            changed = true;
                    }
                }

                if (getFirst(X).count(EPSILON) == 0)
                {
                    allNullable = false;
                    break;
                }
            }

            if (allNullable)
            {
                if (first[A].insert(EPSILON).second)
                    changed = true;
            }
        }
    }
}

void FirstFollowCalculator::computeFollow()
{
    // 起始符加 $
    follow[grammar.getStartSymbol()].insert(END_MARK);

    bool changed = true;
    while (changed)
    {
        changed = false;

        for (const Production &p : grammar.getProductions())
        {
            const Symbol &A = p.left;
            const vector<Symbol> &alpha = p.right;

            for (size_t i = 0; i < alpha.size(); ++i)
            {
                const Symbol &B = alpha[i];
                if (B.type != SymbolType::NON_TERMINAL)
                    continue;

                // β = alpha[i+1 ...]
                vector<Symbol> beta;
                for (size_t j = i + 1; j < alpha.size(); ++j)
                    beta.push_back(alpha[j]);

                set<Symbol> firstBeta = firstOfSequence(beta);

                // FIRST(β) - ε
                for (const Symbol &s : firstBeta)
                {
                    if (s != EPSILON)
                    {
                        if (follow[B].insert(s).second)
                            changed = true;
                    }
                }

                // β ⇒* ε
                if (beta.empty() || firstBeta.count(EPSILON))
                {
                    for (const Symbol &s : follow[A])
                    {
                        if (follow[B].insert(s).second)
                            changed = true;
                    }
                }
            }
        }
    }
}

void Analysis()
{
    string prog;
    read_prog(prog);
    /* 骚年们 请开始你们的表演 */
    /********* Begin *********/

    /********* End *********/
}