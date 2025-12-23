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

const Symbol EPSILON("E", SymbolType::EPSILON);
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

// ========= GrammarBuilder =========

class GrammarBuilder
{
public:
    GrammarBuilder(Grammar &g);

    // 从多行文本读取文法
    void loadFromText(const string &text);

private:
    Grammar &grammar;

    // 符号表，保证同名符号唯一
    map<string, Symbol> symbolTable;

    // 工具函数
    Symbol getOrCreateSymbol(const string &name, SymbolType type);

    void parseLine(const string &line);

    // 字符串工具
    static string trim(const string &s);
    static vector<string> split(const string &s, char delim);
};

// ========= InputProcessor =========

class InputProcessor
{
public:
    InputProcessor(const string &input, const Grammar &grammar);

    // 返回解析用的符号流（末尾含 $）
    vector<Symbol> getSymbolStream() const;

private:
    vector<Symbol> symbolStream;

    // 工具函数
    Symbol makeTerminal(const string &token, const Grammar &grammar);
};

// ========= ParseTreeNode =========

class ParseTreeNode
{
public:
    Symbol symbol;                    // 当前结点对应的文法符号
    vector<ParseTreeNode *> children; // 子结点（从左到右）

    explicit ParseTreeNode(const Symbol &s);

    bool isLeaf() const;
    void addChild(ParseTreeNode *child);
};

// ========= ParseTree =========

class ParseTree
{
public:
    explicit ParseTree(ParseTreeNode *root);

    ParseTreeNode *getRoot() const;
    void print() const;

private:
    ParseTreeNode *root;

    void printNode(ParseTreeNode *node, int indent) const;
};

// ========= LL1Parser =========

class LL1Parser
{
public:
    LL1Parser(const Grammar &g, const LL1Table &table);

    void parse(const std::vector<Symbol> &input);

private:
    const Grammar &grammar;
    const LL1Table &ll1Table;

    std::vector<Symbol> tokens;
    size_t pos;

    ParseTree *parseTree;

    ParseTreeNode *parseSymbol(const Symbol &X);

    void reportError(int line, const std::string &msg);
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

    followSet[END_MARK]; // 保证存在（通常不用，但安全）
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

LL1Parser::LL1Parser(const Grammar &g, const LL1Table &table)
    : grammar(g), ll1Table(table), pos(0), parseTree(nullptr)
{
}

void LL1Parser::parse(const std::vector<Symbol> &input)
{
    tokens = input;
    tokens.push_back(END_MARK);
    pos = 0;

    std::cout << "=== 开始解析 ===\n";

    ParseTreeNode *root = parseSymbol(grammar.getStartSymbol());
    parseTree = new ParseTree(root);

    if (tokens[pos] != END_MARK)
    {
        reportError(0, "输入未完全匹配");
    }

    std::cout << "=== 解析结束 ===\n";
    parseTree->print();
}

ParseTreeNode *LL1Parser::parseSymbol(const Symbol &X)
{
    Symbol a = tokens[pos];
    ParseTreeNode *node = new ParseTreeNode(X);

    // ① 终结符 or $
    if (X.type == SymbolType::TERMINAL || X.type == SymbolType::END)
    {
        if (X == a)
        {
            std::cout << "匹配终结符: " << X.name << " ✓\n";
            pos++;
        }
        else
        {
            reportError(0, "缺少 \"" + X.name + "\"");
        }
        return node;
    }

    // ② ε
    if (X.type == SymbolType::EPSILON)
    {
        return node;
    }

    // ③ 非终结符
    if (X.type == SymbolType::NON_TERMINAL)
    {
        if (!ll1Table.hasEntry(X, a))
        {
            reportError(0, "非法符号 \"" + a.name + "\"");
            return node;
        }

        const Production &p = ll1Table.getEntry(X, a);

        std::cout << "使用产生式: " << p.left.name << " -> ";
        for (const Symbol &s : p.right)
            std::cout << s.name << " ";
        std::cout << std::endl;

        for (const Symbol &Y : p.right)
        {
            ParseTreeNode *child = parseSymbol(Y);
            node->addChild(child);
        }

        return node;
    }

    return node;
}

void LL1Parser::reportError(int line, const std::string &msg)
{
    std::cout << "语法错误, 第" << line << "行, " << msg << std::endl;
}

GrammarBuilder::GrammarBuilder(Grammar &g) : grammar(g)
{
    // 已有预注册 ε 和 $
    symbolTable["ε"] = EPSILON;
    symbolTable["E"] = EPSILON;
    symbolTable["$"] = END_MARK;

    // 预注册保留字和操作符为终结符
    vector<string> terminals = {
        "if", "then", "else", "while",
        "ID", "NUM",
        "{", "}", "(", ")", ";",
        "+", "-", "*", "/",
        "=", "<", "<=", "==", ">", ">="};

    for (const string &t : terminals)
    {
        symbolTable[t] = Symbol(t, SymbolType::TERMINAL);
    }
}

void GrammarBuilder::loadFromText(const string &text)
{
    istringstream iss(text);
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
    // 形如：A -> B c | d
    size_t pos = line.find("->");
    if (pos == string::npos)
        return;

    string leftStr = trim(line.substr(0, pos));
    string rightStr = trim(line.substr(pos + 2));

    // 左部一定是非终结符
    Symbol left = getOrCreateSymbol(leftStr, SymbolType::NON_TERMINAL);

    // 右部按 |
    vector<string> alternatives = split(rightStr, '|');

    for (string &alt : alternatives)
    {
        vector<Symbol> right;

        alt = trim(alt);
        if (alt == "ε" || alt == "E")
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

                // 优先查 symbolTable
                auto it = symbolTable.find(tok);
                if (it != symbolTable.end())
                {
                    right.push_back(it->second);
                }
                else
                {
                    // 没有注册的 token 当作非终结符
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

// ========= InputProcessor 实现 =========

InputProcessor::InputProcessor(const string &input, const Grammar &grammar)
{
    stringstream ss(input);
    string token;

    while (ss >> token)
    {
        // ① 忽略 ε
        if (token == "ε")
            continue;

        // ② 忽略 Ctrl+Z / EOF / 非法控制符
        if (token.size() == 1 && token[0] == 26) // ASCII 0x1A
            continue;

        Symbol s = makeTerminal(token, grammar);
        symbolStream.push_back(s);
    }

    // 末尾加入 $
    symbolStream.push_back(END_MARK);
}

vector<Symbol> InputProcessor::getSymbolStream() const
{
    return symbolStream;
}

Symbol InputProcessor::makeTerminal(const string &token, const Grammar &grammar)
{
    if (token == "ε")
    {
        throw logic_error("EPSILON must not appear in input stream");
    }

    Symbol s(token, SymbolType::TERMINAL);

    if (!grammar.isTerminal(s))
    {
        throw logic_error("Unknown terminal in input: " + token);
    }

    return s;
}

// ========= ParseTreeNode 实现 =========

ParseTreeNode::ParseTreeNode(const Symbol &s) : symbol(s) {}

bool ParseTreeNode::isLeaf() const
{
    return children.empty();
}

void ParseTreeNode::addChild(ParseTreeNode *child)
{
    children.push_back(child);
}

// ========= ParseTree 实现 =========

ParseTree::ParseTree(ParseTreeNode *r) : root(r) {}

ParseTreeNode *ParseTree::getRoot() const
{
    return root;
}

void ParseTree::print() const
{
    printNode(root, 0);
}

void ParseTree::printNode(ParseTreeNode *node, int indent) const
{
    if (!node)
        return;

    for (int i = 0; i < indent; ++i)
        cout << '\t';

    cout << node->symbol.name << endl;

    for (ParseTreeNode *child : node->children)
    {
        printNode(child, indent + 1);
    }
}

void Analysis()
{
    /* ========= 1. 构建文法 ========= */

    Grammar grammar;
    GrammarBuilder builder(grammar);

    string grammarText = R"(
program -> compoundstmt
stmt -> ifstmt | whilestmt | assgstmt | compoundstmt
compoundstmt -> { stmts }
stmts -> stmt stmts | ε
ifstmt -> if ( boolexpr ) then stmt else stmt
whilestmt -> while ( boolexpr ) stmt
assgstmt -> ID = arithexpr ;
boolexpr -> arithexpr boolop arithexpr
boolop -> < | > | <= | >= | ==
arithexpr -> multexpr arithexprprime
arithexprprime -> + multexpr arithexprprime | - multexpr arithexprprime | ε
multexpr -> simpleexpr multexprprime
multexprprime -> * simpleexpr multexprprime | / simpleexpr multexprprime | ε
simpleexpr -> ID | NUM | ( arithexpr )
)";

    builder.loadFromText(grammarText);
    grammar.setStartSymbol(Symbol("program", SymbolType::NON_TERMINAL));

    cout << "Grammar loaded.\n";

    /* ========= 2. FIRST / FOLLOW ========= */

    FirstFollowCalculator ff(grammar);
    ff.computeFirst();
    ff.computeFollow();

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

    /* ========= 3. 构建 LL(1) 表 ========= */

    LL1Table table;
    table.build(grammar, ff);

    if (table.hasConflict())
    {
        cout << "\n⚠ 文法不是 LL(1)，存在冲突，解析终止。\n";
        return;
    }

    cout << "\nLL(1) table built successfully.\n";

    /* ========= 4. 读取输入 ========= */

    string inputText;
    read_prog(inputText);

    InputProcessor processor(inputText, grammar);
    vector<Symbol> input = processor.getSymbolStream();

    /* ========= 5. LL(1) 预测分析 + 构建语法树 ========= */

    LL1Parser parser(grammar, table);

    cout << "\nStart parsing...\n";
    parser.parse(input);

    cout << "\nParsing finished.\n";
}
