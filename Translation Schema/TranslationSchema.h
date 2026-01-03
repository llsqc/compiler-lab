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

    Symbol() : name(""), type(SymbolType::TERMINAL) {};
    Symbol(const string &n, SymbolType t) : name(n), type(t) {}

    bool operator==(const Symbol &other) const { return name == other.name && type == other.type; }

    bool operator!=(const Symbol &other) const { return !(*this == other); }

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
    Grammar(const Symbol &start) : startSymbol(start) { nonTerminals.insert(start); }
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
    bool isTerminal(const Symbol &s) const { return terminals.count(s) > 0; }
    bool isNonTerminal(const Symbol &s) const { return nonTerminals.count(s) > 0; }

    // 获取集合
    const set<Symbol> &getTerminals() const { return terminals; }
    const set<Symbol> &getNonTerminals() const { return nonTerminals; }

private:
    Symbol startSymbol;
    vector<Production> productions;
    map<Symbol, vector<Production>> prodMap;

    set<Symbol> terminals;
    set<Symbol> nonTerminals;
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

// ========= LL1Table =========

class LL1Table
{
public:
    LL1Table() : conflict(false) {}

    void build(const Grammar &g, const FirstFollowCalculator &ff);

    bool hasEntry(const Symbol &nonTerminal, const Symbol &terminal) const { return table.count({nonTerminal, terminal}) > 0; }

    const Production &getEntry(const Symbol &nonTerminal, const Symbol &terminal) const;

    bool hasConflict() const { return conflict; }

private:
    map<pair<Symbol, Symbol>, Production> table;
    bool conflict;
};

// ========= Token =========

struct Token
{
    Symbol symbol; // 终结符类型：ID / INTNUM / REALNUM / +
    string lexeme; // 原始字符串：a / 123 / 3.14
    int line;
};

// ========= InputProcessor =========

class InputProcessor
{
public:
    InputProcessor(const string &input, const Grammar &grammar);

    // 返回解析用的符号流（末尾含 $）
    vector<Token> getTokenStream() const { return tokenStream; }

private:
    vector<Token> tokenStream;

    // 工具函数
    Symbol makeTerminal(const string &token, const Grammar &grammar);
};

// ========= ParseTreeNode =========

class ParseTreeNode
{
public:
    Symbol symbol;                    // 当前结点对应的文法符号
    vector<ParseTreeNode *> children; // 子结点（从左到右）

    explicit ParseTreeNode(const Symbol &s) : symbol(s) {};

    bool isLeaf() const { return children.empty(); };
    void addChild(ParseTreeNode *child) { children.push_back(child); };
};

// ========= ParseTree =========

class ParseTree
{
public:
    explicit ParseTree(ParseTreeNode *root) : root(root) {};

    ParseTreeNode *getRoot() const { return root; };
    void print() const { printNode(root, 0); };

private:
    ParseTreeNode *root;

    void printNode(ParseTreeNode *node, int indent) const;
};

// ========= LL1Parser =========

class LL1Parser
{
public:
    LL1Parser(const Grammar &g, const LL1Table &table) : grammar(g), ll1Table(table), pos(0), parseTree(nullptr) {}

    void parse(const vector<Token> &input);
    ParseTree *getParseTree() const { return parseTree; }

private:
    const Grammar &grammar;
    const LL1Table &ll1Table;

    vector<Token> tokens;
    size_t pos;

    ParseTree *parseTree;

    ParseTreeNode *parseSymbol(const Symbol &X);

    void reportError(int line, const string &msg);
};

enum class ValueType
{
    INT,
    REAL
};

struct ExprResult
{
    ValueType type;
    double value;
    bool hasValue;

    ExprResult() : type(ValueType::INT), value(0), hasValue(false) {}
    ExprResult(ValueType t, double v) : type(t), value(v), hasValue(true) {}
};

struct SymbolEntry
{
    ValueType type;
    double value;
    bool hasValue;

    SymbolEntry() : type(ValueType::INT), value(0), hasValue(false) {}
    SymbolEntry(ValueType t, double v, bool hv) : type(t), value(v), hasValue(hv) {}
};

class SemanticAnalyzer
{
public:
    SemanticAnalyzer(ParseTree *tree, map<string, SymbolEntry> &symtab) : tree(tree), symtab(symtab) {}
    void analyze() { analyzeNode(tree->getRoot()); }

private:
    ParseTree *tree;
    map<string, SymbolEntry> &symtab;

    void analyzeNode(ParseTreeNode *node);

    void handleAssignStmt(ParseTreeNode *node);
    void handleIfStmt(ParseTreeNode *node);
    void handleWhileStmt(ParseTreeNode *node);
    void handleDecls(ParseTreeNode *node);
    void handleDecl(ParseTreeNode *node);

    ExprResult handleArithExpr(ParseTreeNode *node);
    ExprResult handleMultiExpr(ParseTreeNode *node);
    ExprResult handleSimpleExpr(ParseTreeNode *node);

    bool evalBoolExpr(ParseTreeNode *node);

    [[noreturn]] void semanticError(const string &msg);
    static ParseTreeNode *findChild(ParseTreeNode *node, const string &name)
    {
        for (auto c : node->children)
            if (c->symbol.name == name)
                return c;
        return nullptr;
    }

    static vector<ParseTreeNode *> findChildren(ParseTreeNode *node, const string &name)
    {
        vector<ParseTreeNode *> res;
        for (auto c : node->children)
            if (c->symbol.name == name)
                res.push_back(c);
        return res;
    }
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

#pragma region LL1Table
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

const Production &LL1Table::getEntry(const Symbol &nonTerminal, const Symbol &terminal) const
{
    auto it = table.find({nonTerminal, terminal});
    if (it == table.end())
    {
        throw runtime_error("LL(1) table entry not found");
    }
    return it->second;
}
#pragma endregion

#pragma region LL1Parser
void LL1Parser::parse(const vector<Token> &input)
{
    tokens = input;
    pos = 0;
    ParseTreeNode *root = parseSymbol(grammar.getStartSymbol());
    parseTree = new ParseTree(root);
#ifdef DEBUG
    if (tokens[pos].symbol != END_MARK)
    {
        reportError(0, "输入未完全匹配");
    }
    parseTree->print();
#endif
}

ParseTreeNode *LL1Parser::parseSymbol(const Symbol &X)
{
    Symbol a = tokens[pos].symbol;
    int line = tokens[pos].line;
    ParseTreeNode *node = new ParseTreeNode(X);

    // ① 终结符 or $
    if (X.type == SymbolType::TERMINAL || X.type == SymbolType::END)
    {
        if (X == a)
        {
            pos++; // 正常匹配
        }
        else
        {
            // 缺少终结符：报告错误，但不吃输入
            reportError(line, "缺少\"" + X.name + "\"");
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
        // 如果 LL(1) 表中没有入口
        if (!ll1Table.hasEntry(X, a))
        {
            const auto &prods = grammar.getProductionsOf(X);
            for (const Production &p : prods)
            {
                if (p.right.size() == 1 && p.right[0].type == SymbolType::EPSILON)
                {
                    // 选择 ε 产生式
                    ParseTreeNode *eps = new ParseTreeNode(EPSILON);
                    node->addChild(eps);
                    return node;
                }
            }

            reportError(line, "非法符号 \"" + a.name + "\"");
            pos++; // 丢弃一个输入符号
            return node;
        }

        const Production &p = ll1Table.getEntry(X, a);
        for (const Symbol &Y : p.right)
        {
            ParseTreeNode *child = parseSymbol(Y);
            node->addChild(child);
        }

        return node;
    }

    return node;
}

void LL1Parser::reportError(int line, const string &msg)
{
    cout << "语法错误,第" << line << "行," << msg << endl;
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

#pragma region InputProcessor
InputProcessor::InputProcessor(const string &input, const Grammar &grammar)
{
    int line = 1;
    size_t i = 0;

    while (i < input.size())
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

        // 读一个 token（以空格分隔）
        string tok;
        while (i < input.size() && !isspace(input[i]))
        {
            tok += input[i++];
        }

        // ===== 1. 关键字 & 符号 =====
        static set<string> keywords = {
            "int", "real", "if", "then", "else", "while",
            "{", "}", "(", ")", ";",
            "+", "-", "*", "/",
            "=", "<", "<=", "==", ">", ">="};

        if (keywords.count(tok))
        {
            tokenStream.push_back({Symbol(tok, SymbolType::TERMINAL), tok, line});
        }
        // ===== 2. 数字 =====
        else if (isdigit(tok[0]))
        {
            if (tok.find('.') != string::npos)
            {
                tokenStream.push_back({Symbol("REALNUM", SymbolType::TERMINAL), tok, line}); // REALNUM
            }
            else
            {
                tokenStream.push_back({Symbol("INTNUM", SymbolType::TERMINAL), tok, line}); // INTNUM
            }
        }
        // ===== 3. 标识符 =====
        else
        {
            tokenStream.push_back({Symbol("ID", SymbolType::TERMINAL), tok, line}); // ID
        }
    }

    // 结束符 $
    tokenStream.push_back({END_MARK, "$", line});
}

Symbol InputProcessor::makeTerminal(const string &token, const Grammar &grammar)
{
    if (token == "E")
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
#pragma endregion

#pragma region ParseTree
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
#pragma endregion

#pragma region SemanticAnalyzer
void SemanticAnalyzer::analyzeNode(ParseTreeNode *node)
{
    if (!node)
        return;

    const string &name = node->symbol.name;

    if (name == "program")
    {
        // program → decls compoundstmt
        handleDecls(node->children[0]);
        analyzeNode(node->children[1]);
    }
    else if (name == "decls")
    {
        handleDecls(node);
    }
    else if (name == "compoundstmt")
    {
        // { stmts }
        analyzeNode(node->children[1]);
    }
    else if (name == "stmts")
    {
        // stmt stmts | E
        if (node->children.size() == 2)
        {
            analyzeNode(node->children[0]);
            analyzeNode(node->children[1]);
        }
    }
    else if (name == "assgstmt")
    {
        handleAssignStmt(node);
    }
    else if (name == "ifstmt")
    {
        handleIfStmt(node);
    }
    else if (name == "whilestmt")
    {
        handleWhileStmt(node);
    }
}

void SemanticAnalyzer::handleDecls(ParseTreeNode *node)
{
    // decls → decl ; decls | E
    if (node->children.empty())
        return;

    handleDecl(node->children[0]);
    handleDecls(node->children[2]);
}

void SemanticAnalyzer::handleDecl(ParseTreeNode *node)
{
    // decl → int ID = INTNUM
    // decl → real ID = REALNUM

    string typeName = node->children[0]->symbol.name;
    string varName = node->children[1]->children[0]->symbol.name;
    string literal = node->children[3]->children[0]->symbol.name;

    if (symtab.count(varName))
        semanticError("variable redeclared: " + varName);

    if (typeName == "int")
    {
        if (literal.find('.') != string::npos)
            semanticError("realnum can not be translated into int type");

        int v = stoi(literal);
        symtab[varName] = SymbolEntry(ValueType::INT, v, true);
    }
    else if (typeName == "real")
    {
        double v = stod(literal);
        symtab[varName] = SymbolEntry(ValueType::REAL, v, true);
    }
}

void SemanticAnalyzer::handleAssignStmt(ParseTreeNode *node)
{
    string varName = node->children[0]->children[0]->symbol.name;

    if (!symtab.count(varName))
        semanticError("variable not declared: " + varName);

    ExprResult rhs = handleArithExpr(node->children[2]);
    SymbolEntry &entry = symtab[varName];

    if (entry.type == ValueType::INT && rhs.type == ValueType::REAL)
        semanticError("realnum can not be translated into int type");

    entry.value = rhs.value;
    entry.hasValue = true;
}

void SemanticAnalyzer::handleIfStmt(ParseTreeNode *node)
{
    bool cond = evalBoolExpr(node->children[2]);

    if (cond)
        analyzeNode(node->children[5]); // then stmt
    else
        analyzeNode(node->children[7]); // else stmt
}

void SemanticAnalyzer::handleWhileStmt(ParseTreeNode *node)
{
    while (evalBoolExpr(node->children[2]))
    {
        analyzeNode(node->children[4]);
    }
}

bool SemanticAnalyzer::evalBoolExpr(ParseTreeNode *node)
{
    ExprResult lhs = handleArithExpr(node->children[0]);
    ExprResult rhs = handleArithExpr(node->children[2]);

    string op = node->children[1]->children[0]->symbol.name;

    if (op == "<")
        return lhs.value < rhs.value;
    if (op == ">")
        return lhs.value > rhs.value;
    if (op == "<=")
        return lhs.value <= rhs.value;
    if (op == ">=")
        return lhs.value >= rhs.value;
    if (op == "==")
        return lhs.value == rhs.value;

    semanticError("unknown bool operator");
}

ExprResult SemanticAnalyzer::handleArithExpr(ParseTreeNode *node)
{
    ExprResult left = handleMultiExpr(node->children[0]);
    ParseTreeNode *prime = node->children[1];

    while (!prime->children.empty())
    {
        string op = prime->children[0]->symbol.name;
        ExprResult right = handleMultiExpr(prime->children[1]);

        if (op == "+")
            left.value += right.value;
        else
            left.value -= right.value;

        if (left.type == ValueType::INT && right.type == ValueType::REAL)
            left.type = ValueType::REAL;

        prime = prime->children[2];
    }
    return left;
}

ExprResult SemanticAnalyzer::handleMultiExpr(ParseTreeNode *node)
{
    ExprResult left = handleSimpleExpr(node->children[0]);
    ParseTreeNode *prime = node->children[1];

    while (!prime->children.empty())
    {
        string op = prime->children[0]->symbol.name;
        ExprResult right = handleSimpleExpr(prime->children[1]);

        if (op == "*")
            left.value *= right.value;
        else
            left.value /= right.value;

        if (left.type == ValueType::INT && right.type == ValueType::REAL)
            left.type = ValueType::REAL;

        prime = prime->children[2];
    }
    return left;
}

ExprResult SemanticAnalyzer::handleSimpleExpr(ParseTreeNode *node)
{
    ParseTreeNode *child = node->children[0];

    if (child->symbol.name == "ID")
    {
        string name = child->children[0]->symbol.name;

        if (!symtab.count(name))
            semanticError("variable not declared: " + name);

        SymbolEntry &e = symtab[name];
        if (!e.hasValue)
            semanticError("variable not initialized: " + name);

        return ExprResult(e.type, e.value);
    }

    if (child->symbol.name == "INTNUM")
        return ExprResult(ValueType::INT, stod(child->children[0]->symbol.name));

    if (child->symbol.name == "REALNUM")
        return ExprResult(ValueType::REAL, stod(child->children[0]->symbol.name));

    // ( arithexpr )
    return handleArithExpr(child->children[1]);
}

[[noreturn]] void SemanticAnalyzer::semanticError(const string &msg)
{
    cerr << "error message:" << msg << endl;
    exit(1);
}
#pragma endregion

void Analysis()
{
    /* ========= 1. 构建文法 ========= */
    string grammarText = R"(
program -> decls compoundstmt
decls -> decl ; decls | E
decl -> int ID = INTNUM | real ID = REALNUM
stmt -> ifstmt | assgstmt | compoundstmt
compoundstmt -> { stmts }
stmts -> stmt stmts | E
ifstmt -> if ( boolexpr ) then stmt else stmt
assgstmt -> ID = arithexpr ;
boolexpr -> arithexpr boolop arithexpr
boolop -> < | > | <= | >= | ==
arithexpr -> multexpr arithexprprime
arithexprprime -> + multexpr arithexprprime | - multexpr arithexprprime | E
multexpr -> simpleexpr multexprprime
multexprprime -> * simpleexpr multexprprime | / simpleexpr multexprprime | E
simpleexpr -> ID | INTNUM | REALNUM | ( arithexpr )
)";

    vector<string> terminals = {
        "int", "real",
        "if", "then", "else", "while",
        "ID", "INTNUM", "REALNUM",
        "{", "}", "(", ")", ";",
        "+", "-", "*", "/",
        "=", "<", "<=", "==", ">", ">="};

    Grammar grammar;
    GrammarBuilder builder(grammar, terminals, grammarText);

    grammar.setStartSymbol(Symbol("program", SymbolType::NON_TERMINAL));
#ifdef DEBUG
    cout << "Grammar loaded.\n";
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
#endif

    /* ========= 3. 构建 LL(1) 表 ========= */

    LL1Table table;
    table.build(grammar, ff);
#ifdef DEBUG
    if (table.hasConflict())
    {
        cout << "\n⚠ 文法不是 LL(1)，存在冲突，解析终止。\n";
        return;
    }
    cout << "\nLL(1) table built successfully.\n";
#endif

    /* ========= 4. 读取输入 ========= */

    string inputText;
    read_prog(inputText);

    InputProcessor processor(inputText, grammar);
    vector<Token> input = processor.getTokenStream();

    /* ========= 5. LL(1) 预测分析 + 构建语法树 ========= */

    LL1Parser parser(grammar, table);
#ifdef DEBUG
    cout << "\nStart parsing...\n";
#endif
    parser.parse(input);
#ifdef DEBUG
    cout << "\nParsing finished.\n";
#endif

    /* ========= 6. 语义分析 ========= */

    ParseTree *tree = parser.getParseTree();
    map<string, SymbolEntry> symbolTable;
    SemanticAnalyzer semantic(tree, symbolTable);
    semantic.analyze();

    /* ========= 7. 输出最终结果 ========= */
    cout << "\n====== Program Result ======\n";
    for (const auto &p : symbolTable)
    {
        cout << p.first << " = " << p.second.value << endl;
    }
}