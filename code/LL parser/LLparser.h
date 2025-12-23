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

    map<Symbol, set<Symbol>> firstSet;
    map<Symbol, set<Symbol>> followSet;

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

    void synchronize(const Symbol &nonTerminal, const vector<Symbol> &input, size_t &pos);
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
    return !(*this == other);
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
        firstSet[nt] = set<Symbol>();
        followSet[nt] = set<Symbol>();
    }

    // 终结符的 FIRST 是它自己
    for (const Symbol &t : grammar.getTerminals())
    {
        firstSet[t].insert(t);
    }

    firstSet[END_MARK].insert(END_MARK);
}

const set<Symbol> &FirstFollowCalculator::getFirst(const Symbol &s) const
{
    auto it = firstSet.find(s);
    if (it == firstSet.end())
        throw logic_error("FIRST not initialized");
    return it->second;
}

const set<Symbol> &FirstFollowCalculator::getFollow(const Symbol &s) const
{
    auto it = followSet.find(s);
    if (it == followSet.end())
        throw logic_error("Follow not initialized");
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
                        if (firstSet[A].insert(s).second)
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
                if (firstSet[A].insert(EPSILON).second)
                    changed = true;
            }
        }
    }
}

void FirstFollowCalculator::computeFollow()
{
    // 起始符加 $
    followSet[grammar.getStartSymbol()].insert(END_MARK);

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
                        if (followSet[B].insert(s).second)
                            changed = true;
                    }
                }

                // β ⇒* ε
                if (beta.empty() || firstBeta.count(EPSILON))
                {
                    for (const Symbol &s : followSet[A])
                    {
                        if (followSet[B].insert(s).second)
                            changed = true;
                    }
                }
            }
        }
    }
}

// ========= 实现 =========

LL1Table::LL1Table() : conflict(false) {}

void LL1Table::build(const Grammar &g, const FirstFollowCalculator &ff)
{
    conflict = false;
    table.clear();

    for (const Production &p : g.getProductions())
    {
        const Symbol &A = p.left;
        const vector<Symbol> &alpha = p.right;

        // FIRST(α)
        set<Symbol> firstAlpha = ff.firstOfSequence(alpha);

        // 情况 1：a ∈ FIRST(α) - {ε}
        for (const Symbol &a : firstAlpha)
        {
            if (a.type == SymbolType::EPSILON)
                continue;

            pair<Symbol, Symbol> key = {A, a};

            if (table.count(key))
            {
                conflict = true;
            }
            else
            {
                table[key] = p;
            }
        }

        // 情况 2：ε ∈ FIRST(α)
        if (firstAlpha.count(EPSILON))
        {
            const set<Symbol> &followA = ff.getFollow(A);

            for (const Symbol &b : followA)
            {
                pair<Symbol, Symbol> key = {A, b};

                if (table.count(key))
                {
                    conflict = true;
                }
                else
                {
                    table[key] = p;
                }
            }
        }
    }
}

bool LL1Table::hasEntry(const Symbol &nonTerminal, const Symbol &terminal) const
{
    return table.count({nonTerminal, terminal}) > 0;
}

const Production &LL1Table::getEntry(const Symbol &nonTerminal, const Symbol &terminal) const
{
    auto it = table.find({nonTerminal, terminal});
    if (it == table.end())
    {
        throw runtime_error("LL(1) table entry not found");
    }
    return it->second;
}

bool LL1Table::hasConflict() const
{
    return conflict;
}

LL1Parser::LL1Parser(const Grammar &g, const LL1Table &table) : grammar(g), ll1Table(table)
{
}

void LL1Parser::reportError(int line, const string &msg)
{
    cout << "语法错误,第" << line << "行," << msg << endl;
}

void LL1Parser::synchronize(const Symbol &nonTerminal, const vector<Symbol> &input, size_t &pos)
{
    // 跳过输入符号，直到：
    // 1) 遇到 FOLLOW(nonTerminal) 中的符号
    // 2) 或遇到输入结束

    // 这里假设 FOLLOW 集已经在表构造阶段间接使用过
    while (pos < input.size())
    {
        // 如果表中存在 (nonTerminal, 当前输入)，可以继续分析
        if (ll1Table.hasEntry(nonTerminal, input[pos]))
            return;

        pos++; // 丢弃一个输入符号
    }
}

void LL1Parser::parse(const vector<Symbol> &input)
{
    parseStack = stack<Symbol>();

    // 初始化分析栈
    parseStack.push(END_MARK);
    parseStack.push(grammar.getStartSymbol());

    size_t pos = 0;

    // 输入末尾补 $
    vector<Symbol> tokens = input;
    tokens.push_back(END_MARK);

    while (!parseStack.empty())
    {
        Symbol X = parseStack.top();
        Symbol a = tokens[pos];

        // 1️⃣ 栈顶是终结符 或 $
        if (X.type == SymbolType::TERMINAL || X.type == SymbolType::END)
        {
            if (X == a)
            {
                parseStack.pop();
                pos++;
            }
            else
            {
                reportError(0, "缺少 \"" + X.name + "\"");
                parseStack.pop(); // 弹出栈顶，试图继续
            }
        }
        // 2️⃣ 栈顶是 ε
        else if (X.type == SymbolType::EPSILON)
        {
            parseStack.pop(); // ε 直接弹出
        }
        // 3️⃣ 栈顶是非终结符
        else if (X.type == SymbolType::NON_TERMINAL)
        {
            if (ll1Table.hasEntry(X, a))
            {
                const Production &p = ll1Table.getEntry(X, a);
                parseStack.pop();

                // 逆序压入产生式右部
                for (auto it = p.right.rbegin(); it != p.right.rend(); ++it)
                {
                    if (it->type != SymbolType::EPSILON)
                        parseStack.push(*it);
                }
            }
            else
            {
                reportError(0, "非法符号 \"" + a.name + "\"");
                parseStack.pop();
                synchronize(X, tokens, pos);
            }
        }
        else
        {
            // 理论上不会到这里
            parseStack.pop();
        }
    }

    if (pos < tokens.size() - 1)
    {
        reportError(0, "输入未完全匹配");
    }
}

void Analysis()
{
    /********* Begin *********/

    // ---------- 1. 定义符号 ----------
    Symbol Program("Program", SymbolType::NON_TERMINAL);
    Symbol Stmt("Stmt", SymbolType::NON_TERMINAL);

    Symbol ID("ID", SymbolType::TERMINAL);
    Symbol ASSIGN("=", SymbolType::TERMINAL);
    Symbol NUM("NUM", SymbolType::TERMINAL);
    Symbol SEMI(";", SymbolType::TERMINAL);

    // ---------- 2. 构建文法 ----------
    Grammar grammar(Program);

    // Program → Stmt
    grammar.addProduction(Program, {Stmt});

    // Stmt → ID = NUM ;
    grammar.addProduction(Stmt, {ID, ASSIGN, NUM, SEMI});

    // ---------- 3. FIRST / FOLLOW ----------
    FirstFollowCalculator ff(grammar);
    ff.computeFirst();
    ff.computeFollow();

    // ---------- 4. 构建 LL(1) 表 ----------
    LL1Table table;
    table.build(grammar, ff);

    if (table.hasConflict())
    {
        cout << "该文法不是 LL(1) 文法" << endl;
        return;
    }

    // ---------- 5. 构造输入 ----------
    vector<Symbol> input = {
        ID, ASSIGN, NUM, SEMI};

    // ---------- 6. 预测分析 ----------
    LL1Parser parser(grammar, table);
    parser.parse(input);

    cout << "分析成功" << endl;

    /********* End *********/
}
