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
#include <algorithm>
using namespace std;
// #define DEBUG
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

// ========= Symbol =========

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

    Symbol() : name(""), type(SymbolType::TERMINAL) {};
    Symbol(const string &n, SymbolType t) : name(n), type(t) {}

    bool operator==(const Symbol &other) const
    {
        return name == other.name && type == other.type;
    }

    bool operator!=(const Symbol &other) const
    {
        return !(*this == other);
    }

    bool operator<(const Symbol &other) const
    {
        if (type != other.type)
            return type < other.type;
        return name < other.name;
    }
};

// ========= 规范ε和$符号 =========

const Symbol EPSILON("E", SymbolType::EPSILON);
const Symbol END_MARK("$", SymbolType::END);

// ========= Production =========

struct Production
{
    Symbol left;          // A
    vector<Symbol> right; // α

    // A → α
    Production() {}
    Production(const Symbol &l, const vector<Symbol> &r) : left(l), right(r) {}
};

// ========= Grammar =========

class Grammar
{
public:
    // 构造
    Grammar() {};
    Grammar(const Symbol &start) { nonTerminals.insert(start); }

    void setStartSymbol(const Symbol &s)
    {
        startSymbol = s;
        nonTerminals.insert(s);
    }
    Symbol getStartSymbol() const { return startSymbol; }

    // 添加产生式 A → α
    void addProduction(const Symbol &left, const vector<Symbol> &right);

    // 查询
    const vector<Production> &getProductions() const { return productions; }
    const vector<Production> &getProductionsOf(const Symbol &nt) const;

    // 判断符号类型
    bool isTerminal(const Symbol &s) const { return terminals.find(s) != terminals.end(); }
    bool isNonTerminal(const Symbol &s) const { return nonTerminals.find(s) != nonTerminals.end(); }

    // 获取集合
    const set<Symbol> &getTerminals() const { return terminals; }
    const set<Symbol> &getNonTerminals() const { return nonTerminals; }

    // 添加和获取增广文法
    void augmentGrammar();
    const Symbol &getAugmentedStart() const { return augmentedStartSymbol; }
    int getAugmentedProductionIndex() const { return augmentedProductionIndex; }

    // 获取符号，确保同名唯一
    const Symbol &getTerminalByName(const string &name) const;
    const Symbol &getNonTerminalByName(const string &name) const;

private:
    Symbol startSymbol;
    vector<Production> productions;
    map<Symbol, vector<Production>> prodMap;

    set<Symbol> terminals;
    set<Symbol> nonTerminals;

    Symbol augmentedStartSymbol;
    int augmentedProductionIndex = -1;
};

// ========= GrammarBuilder =========

class GrammarBuilder
{
public:
    GrammarBuilder(Grammar &g, const vector<string> &terminals, const string &grammarText);

private:
    Grammar &grammar;

    // 符号表，保证同名唯一
    map<string, Symbol> symbolTable;

    // 内部工具
    Symbol getOrCreateSymbol(const string &name, SymbolType type);
    void parseLine(const string &line);

    // 字符串工具
    static string trim(const string &s);
    static vector<string> split(const string &s, char delim);
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
};

// ========= Token =========

struct Token
{
    Symbol symbol;
    int line;
};

// ========= InputProcessor =========

class InputProcessor
{
public:
    InputProcessor(const string &input, const Grammar &grammar);

    // 返回解析用的符号流（末尾含 $）
    vector<Token> getTokenStream() const;

private:
    vector<Token> tokenStream;

    // 工具函数
    Symbol makeTerminal(const string &token, const Grammar &grammar);
};

// ========= LRItem =========

class LRItem
{
public:
    // 指向 Grammar::productions 中的产生式下标
    int productionIndex;

    // 点的位置：0 <= dotPos <= right.size()
    int dotPos;

    LRItem(int prodIdx, int dot) : productionIndex(prodIdx), dotPos(dot) {}

    bool operator==(const LRItem &other) const
    {
        return productionIndex == other.productionIndex && dotPos == other.dotPos;
    }

    bool operator<(const LRItem &other) const
    {
        if (productionIndex != other.productionIndex)
            return productionIndex < other.productionIndex;
        return dotPos < other.dotPos;
    }
};

// ========= ItemSet =========

class ItemSet
{
public:
    set<LRItem> items;

    ItemSet() = default;

    explicit ItemSet(const set<LRItem> &its) : items(its) {}

    bool operator==(const ItemSet &other) const
    {
        return items == other.items;
    }

    bool operator<(const ItemSet &other) const
    {
        return items < other.items;
    }
};

// ========= LR0Automaton =========

class LR0Automaton
{
public:
    explicit LR0Automaton(const Grammar &g) : grammar(g) {}

    void build(); // 构造整个项目集族

    const vector<ItemSet> &getStates() const;
    int getGoto(int state, const Symbol &X) const;

private:
    const Grammar &grammar;

    vector<ItemSet> states; // 状态编号 = 下标
    map<pair<int, Symbol>, int> gotoTable;

    // —— 核心算法 ——
    ItemSet closure(const ItemSet &I);
    ItemSet gotoState(const ItemSet &I, const Symbol &X);

    int findState(const ItemSet &S) const
    {
        for (int i = 0; i < (int)states.size(); ++i)
        {
            if (states[i] == S)
                return i;
        }
        return -1;
    }
};

// ========= Action =========

enum class ActionType
{
    SHIFT,
    REDUCE,
    ACCEPT,
    ERROR
};

struct Action
{
    ActionType type;
    int target; // shift: 状态号 | reduce: 产生式编号

    Action() : type(ActionType::ERROR), target(-1) {}
    Action(ActionType t, int v) : type(t), target(v) {}
};

// ========= SLRParsingTable =========

class SLRParsingTable
{
public:
    SLRParsingTable() : conflict(false) {}

    void build(const Grammar &g, const LR0Automaton &automaton, const FirstFollowCalculator &ff);

    const Action &getAction(int state, const Symbol &terminal) const;
    int getGoto(int state, const Symbol &nonTerminal) const;

    bool hasConflict() const { return conflict; }

private:
    // ACTION[state, terminal]
    map<pair<int, Symbol>, Action> actionTable;

    // GOTO[state, nonterminal]
    map<pair<int, Symbol>, int> gotoTable;

    bool conflict;

    void addAction(int state, const Symbol &a, const Action &act);
};

// ======== SLRParser =========

class SLRParser
{
public:
    SLRParser(const Grammar &g, const SLRParsingTable &table) : grammar(g), parsingTable(table), pos(0) {}

    void parse(const vector<Token> &input);

private:
    const Grammar &grammar;
    const SLRParsingTable &parsingTable;

    vector<Token> tokens;
    size_t pos;

    stack<int> stateStack;
    stack<Symbol> symbolStack;

    vector<int> reduceSequence;

    void reportError(int line, const string &msg);
    void printRightmostDerivation() const;
    string sentenceToString(const vector<Symbol> &sentence) const;
};

#pragma region Grammar
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

const vector<Production> &Grammar::getProductionsOf(const Symbol &nt) const
{
    static vector<Production> empty;
    auto it = prodMap.find(nt);
    if (it == prodMap.end())
        return empty;
    return it->second;
}

void Grammar::augmentGrammar()
{
    // 原始开始符号，例如 program
    Symbol originalStart = startSymbol;

    // 新开始符号 program'
    augmentedStartSymbol =
        Symbol(originalStart.name + "'", SymbolType::NON_TERMINAL);

    // S' -> S
    Production p;
    p.left = augmentedStartSymbol;
    p.right.push_back(originalStart);

    // 插到 productions 的最前面（非常重要）
    productions.insert(productions.begin(), p);

    augmentedProductionIndex = 0;

    // 更新非终结符集合
    nonTerminals.insert(augmentedStartSymbol);

    // 更新 startSymbol
    startSymbol = augmentedStartSymbol;
}

const Symbol &Grammar::getTerminalByName(const string &name) const
{
    for (const Symbol &t : terminals)
    {
        if (t.name == name)
            return t;
    }
    throw logic_error("Unknown terminal: " + name);
}

const Symbol &Grammar::getNonTerminalByName(const string &name) const
{
    for (const Symbol &nt : nonTerminals)
    {
        if (nt.name == name)
            return nt;
    }
    throw logic_error("Unknown non-terminal: " + name);
}
#pragma endregion

#pragma region GrammarBuilder
GrammarBuilder::GrammarBuilder(Grammar &g, const vector<string> &terminals, const string &grammarText) : grammar(g)
{
    /* ===== 1. 预注册 ε 和 $ ===== */
    symbolTable["E"] = EPSILON;
    symbolTable["$"] = END_MARK;

    /* ===== 2. 注册终结符 ===== */
    for (const string &t : terminals)
    {
        symbolTable[t] = Symbol(t, SymbolType::TERMINAL);
    }

    /* ===== 3. 解析文法文本 ===== */
    istringstream iss(grammarText);
    string line;

    while (getline(iss, line))
    {
        line = trim(line);
        if (line.empty())
            continue;
        parseLine(line);
    }
}

void GrammarBuilder::parseLine(const string &line)
{
    // A -> B c | d
    size_t pos = line.find("->");
    if (pos == string::npos)
        return;

    string leftStr = trim(line.substr(0, pos));
    string rightStr = trim(line.substr(pos + 2));

    // 左部一定是非终结符
    Symbol left = getOrCreateSymbol(leftStr, SymbolType::NON_TERMINAL);

    vector<string> alternatives = split(rightStr, '|');

    for (string &alt : alternatives)
    {
        alt = trim(alt);
        vector<Symbol> right;

        if (alt == "E")
        {
            right.push_back(EPSILON);
        }
        else
        {
            vector<string> tokens = split(alt, ' ');
            for (string &tok : tokens)
            {
                tok = trim(tok);
                if (tok.empty())
                    continue;

                auto it = symbolTable.find(tok);
                if (it != symbolTable.end())
                {
                    right.push_back(it->second);
                }
                else
                {
                    // 默认：未注册的符号视为非终结符
                    right.push_back(getOrCreateSymbol(tok, SymbolType::NON_TERMINAL));
                }
            }
        }

        grammar.addProduction(left, right);
    }
}

Symbol GrammarBuilder::getOrCreateSymbol(const string &name, SymbolType type)
{
    auto it = symbolTable.find(name);
    if (it != symbolTable.end())
        return it->second;

    Symbol s(name, type);
    symbolTable[name] = s;
    return s;
}

string GrammarBuilder::trim(const string &s)
{
    size_t l = s.find_first_not_of(" \t\r\n");
    size_t r = s.find_last_not_of(" \t\r\n");
    if (l == string::npos)
        return "";
    return s.substr(l, r - l + 1);
}

vector<string> GrammarBuilder::split(const string &s, char delim)
{
    vector<string> res;
    string tmp;
    istringstream iss(s);

    while (getline(iss, tmp, delim))
        res.push_back(tmp);

    return res;
}
#pragma endregion

#pragma region FirstFollowCalculator

FirstFollowCalculator::FirstFollowCalculator(const Grammar &g) : grammar(g)
{
    // 非终结符
    for (const Symbol &nt : grammar.getNonTerminals())
    {
        firstSet[nt] = set<Symbol>();
        followSet[nt] = set<Symbol>();
    }

    // 终结符
    for (const Symbol &t : grammar.getTerminals())
    {
        firstSet[t].insert(t);
    }

    // ===== 关键补充 =====
    firstSet[EPSILON].insert(EPSILON);
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
#pragma endregion

#pragma region InputProcessor
InputProcessor::InputProcessor(const string &input, const Grammar &grammar)
{
    string token;
    int line = 0;

    for (size_t i = 0; i < input.size();)
    {
        if (input[i] == '\n')
        {
            line++;
            i++;
            continue;
        }

        if (isspace(input[i]))
        {
            i++;
            continue;
        }

        // 读一个 token
        string tok;
        while (i < input.size() && !isspace(input[i]))
        {
            tok += input[i++];
        }

        if (tok == "E")
            continue;

        Symbol s = makeTerminal(tok, grammar);
        tokenStream.push_back({s, line});
    }

    tokenStream.push_back({END_MARK, line});
}

vector<Token> InputProcessor::getTokenStream() const
{
    return tokenStream;
}

Symbol InputProcessor::makeTerminal(const string &token, const Grammar &grammar)
{
    if (token == "E")
    {
        throw logic_error("EPSILON must not appear in input stream");
    }

    const Symbol &s = grammar.getTerminalByName(token);
    return s;
}

#pragma endregion

#pragma region LR0Automaton
ItemSet LR0Automaton::closure(const ItemSet &I)
{
    ItemSet result = I;
    bool changed = true;

    const auto &productions = grammar.getProductions();

    while (changed)
    {
        changed = false;
        set<LRItem> newItems = result.items;

        for (const LRItem &item : result.items)
        {
            const Production &prod = productions[item.productionIndex];

            if ((size_t)item.dotPos >= prod.right.size())
                continue; // 点在末尾无法扩展

            Symbol B = prod.right[item.dotPos];

            if (!grammar.isNonTerminal(B))
                continue; // 只对非终结符做 closure

            // 对每个 B -> γ
            for (size_t i = 0; i < productions.size(); ++i)
            {
                if (productions[i].left == B)
                {
                    const Production &p = productions[i];

                    // 如果是空产生式，dot 直接在末尾
                    int dotPos = (p.right.size() == 1 && p.right[0] == EPSILON) ? 1 : 0;

                    LRItem newItem{(int)i, dotPos};
                    if (newItems.insert(newItem).second)
                    {
                        changed = true;
                    }
                }
            }
        }
        result.items = newItems;
    }

    return result;
}

ItemSet LR0Automaton::gotoState(const ItemSet &I, const Symbol &X)
{
    if (X.type == SymbolType::EPSILON)
        return ItemSet();

    ItemSet J;

    for (const LRItem &item : I.items)
    {
        const Production &p = grammar.getProductions()[item.productionIndex];
        if (item.dotPos < (int)p.right.size() && p.right[item.dotPos] == X)
        {
            LRItem nextItem{item.productionIndex, item.dotPos + 1};
            J.items.insert(nextItem);
        }
    }

    if (!J.items.empty())
    {
        J = closure(J);
    }

    return J;
}

void LR0Automaton::build()
{
    states.clear();
    gotoTable.clear();

    // 初始状态 S' -> · S
    set<LRItem> startItems;
    startItems.insert(LRItem{0, 0});
    ItemSet I0 = closure(ItemSet{startItems});
    states.push_back(I0);

    for (int i = 0; i < (int)states.size(); ++i)
    {
        const ItemSet I = states[i];
        set<Symbol> symbols;
        symbols.insert(grammar.getNonTerminals().begin(), grammar.getNonTerminals().end());
        symbols.insert(grammar.getTerminals().begin(), grammar.getTerminals().end());

        for (const Symbol &X : symbols)
        {
            if (X.type == SymbolType::EPSILON)
                continue; // ⚠ 忽略 EPSILON
            ItemSet J = gotoState(I, X);
            if (J.items.empty())
                continue;
            int j = findState(J);
            if (j == -1)
            {
                j = states.size();
                states.push_back(J);
            }
            gotoTable[{i, X}] = j;
        }
    }
}

const vector<ItemSet> &LR0Automaton::getStates() const
{
    return states;
}

int LR0Automaton::getGoto(int state, const Symbol &X) const
{
    auto it = gotoTable.find({state, X});
    if (it != gotoTable.end())
    {
        return it->second; // 返回转移后的状态编号
    }
    return -1; // 如果没有转移，返回 -1
}
#pragma endregion

#pragma region SLRParsingTable
void SLRParsingTable::build(const Grammar &g, const LR0Automaton &automaton, const FirstFollowCalculator &ff)
{
    actionTable.clear();
    gotoTable.clear();
    conflict = false;

    const auto &states = automaton.getStates();
    const auto &productions = g.getProductions();

    /* =====================================================
     * 1. 填 ACTION 表中的 SHIFT（完全基于 LR(0) GOTO）
     *    若 goto(i, a) = j 且 a 是终结符
     *    ⇒ ACTION[i, a] = shift j
     * ===================================================== */
    for (int i = 0; i < (int)states.size(); ++i)
    {
        for (const Symbol &a : g.getTerminals())
        {
            int j = automaton.getGoto(i, a);
            if (j != -1)
            {
                addAction(i, a, Action(ActionType::SHIFT, j));
            }
        }
    }

    /* =====================================================
     * 2. 填 ACTION 表中的 REDUCE / ACCEPT
     * ===================================================== */
    for (int i = 0; i < (int)states.size(); ++i)
    {
        const ItemSet &I = states[i];

        for (const LRItem &item : I.items)
        {
            const Production &p = productions[item.productionIndex];

            // ---------- 点在末尾 ----------
            if (item.dotPos == (int)p.right.size())
            {
                // S' → S ·  ⇒ ACCEPT
                if (item.productionIndex == g.getAugmentedProductionIndex())
                {
                    addAction(i, END_MARK, Action(ActionType::ACCEPT, -1));
                }
                // A → α ·  ⇒ 对 FOLLOW(A) 中的符号做 REDUCE
                else
                {
                    const Symbol &A = p.left;
                    const auto &followA = ff.getFollow(A);

                    for (const Symbol &a : followA)
                    {
                        addAction(i, a, Action(ActionType::REDUCE, item.productionIndex));
                    }
                }
            }
        }
    }

    /* =====================================================
     * 3. 填 GOTO 表（非终结符）
     * ===================================================== */
    for (int i = 0; i < (int)states.size(); ++i)
    {
        for (const Symbol &A : g.getNonTerminals())
        {
            int j = automaton.getGoto(i, A);
            if (j != -1)
            {
                gotoTable[{i, A}] = j;
            }
        }
    }
#ifdef DEBUG
    cout << "\n===== ACTION Table =====\n";
    for (auto &[key, act] : actionTable)
    {
        cout << "ACTION[" << key.first << ", " << key.second.name << "] = ";
        if (act.type == ActionType::SHIFT)
            cout << "shift " << act.target;
        else if (act.type == ActionType::REDUCE)
            cout << "reduce " << act.target;
        else if (act.type == ActionType::ACCEPT)
            cout << "accept";
        cout << "\n";
    }
    cout << "\n===== GOTO Table =====\n";
    for (auto &[key, val] : gotoTable)
    {
        cout << "GOTO[" << key.first << ", " << key.second.name
             << "] = " << val << "\n";
    }
#endif
}

void SLRParsingTable::addAction(int state, const Symbol &a, const Action &act)
{
    auto key = make_pair(state, a);

    if (actionTable.count(key))
    {
        conflict = true;
        // 可选：输出冲突信息
        // cerr << "SLR 冲突: state " << state << ", symbol " << a.name << endl;
        return;
    }
    actionTable[key] = act;
}
const Action &SLRParsingTable::getAction(int state, const Symbol &terminal) const
{
    static Action errorAction;
    auto it = actionTable.find({state, terminal});
    if (it != actionTable.end())
        return it->second;
    return errorAction;
}

int SLRParsingTable::getGoto(int state, const Symbol &nonTerminal) const
{
    auto it = gotoTable.find({state, nonTerminal});
    if (it != gotoTable.end())
        return it->second;
    return -1;
}

#pragma endregion

#pragma region SLRParser
void SLRParser::parse(const vector<Token> &input)
{
    // 初始化
    tokens = input;
    pos = 0;

    while (!stateStack.empty())
        stateStack.pop();
    while (!symbolStack.empty())
        symbolStack.pop();
    reduceSequence.clear();

    stateStack.push(0);
    symbolStack.push(END_MARK);

    while (true)
    {
        int state = stateStack.top();
        const Symbol &a = tokens[pos].symbol;
        const Action &action = parsingTable.getAction(state, a);
#ifdef DEBUG
        cout << "\n=== PARSER STEP ===\n";

        // 状态栈
        cout << "State stack: ";
        stack<int> ss = stateStack;
        vector<int> ssv;
        while (!ss.empty())
        {
            ssv.push_back(ss.top());
            ss.pop();
        }
        reverse(ssv.begin(), ssv.end());
        for (int x : ssv)
            cout << x << " ";
        cout << "\n";

        // 符号栈
        cout << "Symbol stack: ";
        stack<Symbol> sy = symbolStack;
        vector<Symbol> syv;
        while (!sy.empty())
        {
            syv.push_back(sy.top());
            sy.pop();
        }
        reverse(syv.begin(), syv.end());
        for (auto &s : syv)
            cout << s.name << " ";
        cout << "\n";

        // 向前看
        cout << "Lookahead: " << a.name << "\n";

        cout << "Action: ";
        if (action.type == ActionType::SHIFT)
            cout << "SHIFT " << action.target << "\n";
        else if (action.type == ActionType::REDUCE)
            cout << "REDUCE " << action.target << "\n";
        else if (action.type == ActionType::ACCEPT)
            cout << "ACCEPT\n";
        else
            cout << "ERROR\n";
#endif
        // ===== SHIFT =====
        if (action.type == ActionType::SHIFT)
        {
            symbolStack.push(a);
            stateStack.push(action.target);
            pos++;
        }
        // ===== REDUCE =====
        else if (action.type == ActionType::REDUCE)
        {
            int prodIdx = action.target;
            const Production &p = grammar.getProductions()[prodIdx];

            // 空产生式不弹栈
            int popCount = (p.right.size() == 1 && p.right[0] == EPSILON) ? 0 : (int)p.right.size();

            for (int i = 0; i < popCount; ++i)
            {
                symbolStack.pop();
                stateStack.pop();
            }

            int topState = stateStack.top();
            int gotoState = parsingTable.getGoto(topState, p.left);

            if (gotoState == -1)
            {
                cerr << "ERROR: Invalid GOTO after REDUCE for production " << prodIdx << endl;
                exit(1);
            }

            symbolStack.push(p.left);
            stateStack.push(gotoState);

            reduceSequence.push_back(prodIdx);
        }
        // ===== ACCEPT =====
        else if (action.type == ActionType::ACCEPT)
        {
            printRightmostDerivation();
            return;
        }
        // ===== ERROR（暂不处理）=====
        else
        {
            reportError(tokens[pos].line, "SLR 错误");
            // 实验此阶段不考虑错误处理，直接退出
            return;
        }
    }
}

void SLRParser::reportError(int line, const string &msg)
{
    cout << "语法错误,第" << line << "行," << msg << endl;
}

string SLRParser::sentenceToString(const vector<Symbol> &sentence) const
{
    string res;
    for (const Symbol &s : sentence)
    {
        res += s.name;
        res += " ";
    }
    return res;
}

void SLRParser::printRightmostDerivation() const
{
    // 使用原始开始符号初始化
    vector<Symbol> sent;
    // 初始化为增广产生式右部，也就是原始开始符号
    const Production &augProd = grammar.getProductions()[grammar.getAugmentedProductionIndex()];
    sent = augProd.right;

    string originalStartName = sent[0].name;
    cout << originalStartName << " =>\n";

    // 逆序使用规约
    for (int i = (int)reduceSequence.size() - 1; i >= 0; --i)
    {
        // 跳过增广文法规约
        if (reduceSequence[i] == grammar.getAugmentedProductionIndex())
            continue;

        const Production &p = grammar.getProductions()[reduceSequence[i]];

        // 从右往左找 p.left
        for (int k = (int)sent.size() - 1; k >= 0; --k)
        {
            if (sent[k] == p.left)
            {
                // 删除 A
                sent.erase(sent.begin() + k);

                // 插入 α（若不是 ε）
                if (!(p.right.size() == 1 && p.right[0] == EPSILON))
                {
                    sent.insert(sent.begin() + k, p.right.begin(), p.right.end());
                }
                break;
            }
        }

        if (i > 0)
        {
            cout << sentenceToString(sent) << " =>\n";
        }
    }
    cout << sentenceToString(sent);
}

#pragma endregion

void Analysis()
{
    /* ========= 1. 构建文法 ========= */
    string grammarText = R"(
program -> compoundstmt
stmt -> ifstmt | whilestmt | assgstmt | compoundstmt
compoundstmt -> { stmts }
stmts -> stmt stmts | E
ifstmt -> if ( boolexpr ) then stmt else stmt
whilestmt -> while ( boolexpr ) stmt
assgstmt -> ID = arithexpr ;
boolexpr -> arithexpr boolop arithexpr
boolop -> < | > | <= | >= | ==
arithexpr -> multexpr arithexprprime
arithexprprime -> + multexpr arithexprprime | - multexpr arithexprprime | E
multexpr -> simpleexpr multexprprime
multexprprime -> * simpleexpr multexprprime | / simpleexpr multexprprime | E
simpleexpr -> ID | NUM | ( arithexpr )
)";

    vector<string> terminals = {
        "if", "then", "else", "while",
        "ID", "NUM",
        "{", "}", "(", ")", ";",
        "+", "-", "*", "/",
        "=", "<", "<=", "==", ">", ">="};

    Grammar grammar;
    GrammarBuilder builder(grammar, terminals, grammarText);

    grammar.setStartSymbol(Symbol("program", SymbolType::NON_TERMINAL));
    grammar.augmentGrammar();
#ifdef DEBUG
    cout << "\n===== Grammar =====\n";
    cout << "Start symbol: " << grammar.getStartSymbol().name << "\n";

    const auto &prods = grammar.getProductions();
    for (size_t i = 0; i < prods.size(); i++)
    {
        cout << i << ": " << prods[i].left.name << " -> ";
        for (auto &s : prods[i].right)
            cout << s.name << " ";
        cout << "\n";
    }
#endif

    /* ========= 2. FIRST / FOLLOW ========= */

    FirstFollowCalculator ff(grammar);
    ff.computeFirst();
    ff.computeFollow();
#ifdef DEBUG
    cout << "\nFIRST 集合:\n";
    for (const Symbol &nt : grammar.getNonTerminals())
    {
        cout << "FIRST(" << nt.name << ") = { ";
        for (const Symbol &s : ff.getFirst(nt))
            cout << s.name << " ";
        cout << "}\n";
    }

    cout << "\nFOLLOW 集合:\n";
    for (const Symbol &nt : grammar.getNonTerminals())
    {
        cout << "FOLLOW(" << nt.name << ") = { ";
        for (const Symbol &s : ff.getFollow(nt))
            cout << s.name << " ";
        cout << "}\n";
    }
    if (!ff.getFollow(grammar.getStartSymbol()).count(END_MARK))
    {
        cout << "⚠ FOLLOW(startSymbol) does NOT contain $\n";
    }
#endif

    /* ========= 3. 读取输入 ========= */

    string inputText;
    read_prog(inputText);

    InputProcessor processor(inputText, grammar);
    vector<Token> input = processor.getTokenStream();
#ifdef DEBUG
    cout << "\n===== Input Tokens =====\n";
    for (auto &t : input)
    {
        cout << t.symbol.name << " (line " << t.line << ")\n";
    }
#endif
    /* ========= 4. 构造 LR(0) 自动机 ========= */

    LR0Automaton automaton(grammar);
    automaton.build();

#ifdef DEBUG
    cout << "\nLR(0) 自动机构造完成，共 " << automaton.getStates().size() << " 个状态\n";
    cout << "\n===== LR(0) States =====\n";
    const auto &states = automaton.getStates();
    for (size_t i = 0; i < states.size(); i++)
    {
        cout << "I" << i << ":\n";
        for (auto &item : states[i].items)
        {
            const auto &p = grammar.getProductions()[item.productionIndex];
            cout << "  " << p.left.name << " -> ";
            for (size_t k = 0; k < p.right.size(); k++)
            {
                if (k == (size_t)item.dotPos)
                    cout << "· ";
                cout << p.right[k].name << " ";
            }
            if ((size_t)item.dotPos == p.right.size())
                cout << "·";
            cout << "\n";
        }
        cout << "\n";
    }

    cout << "\n===== LR(0) GOTO =====\n";
    for (size_t i = 0; i < states.size(); i++)
    {
        for (auto &s : grammar.getTerminals())
        {
            int j = automaton.getGoto(i, s);
            if (j != -1)
                cout << "I" << i << " --" << s.name << "--> I" << j << "\n";
        }
        for (auto &s : grammar.getNonTerminals())
        {
            int j = automaton.getGoto(i, s);
            if (j != -1)
                cout << "I" << i << " --" << s.name << "--> I" << j << "\n";
        }
    }

#endif

    /* ========= 5. 构造 SLR 分析表 ========= */

    SLRParsingTable slrTable;
    slrTable.build(grammar, automaton, ff);

#ifdef DEBUG
    if (slrTable.hasConflict())
        cout << "⚠ SLR 表存在冲突\n";
    else
        cout << "SLR 表构造完成，无冲突\n";
#endif

    /* ========= 6. SLR 语法分析 ========= */

    SLRParser parser(grammar, slrTable);
    parser.parse(input);
}