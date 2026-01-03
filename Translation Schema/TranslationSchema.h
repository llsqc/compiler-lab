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
    InputProcessor(const string &input);

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
    string lexeme;                    // 词素（仅终结符有效）
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

    void analyze()
    {
#ifdef DEBUG
        cout << "\n====== 开始语义分析 ======\n";
#endif
        analyzeNode(tree->getRoot());
#ifdef DEBUG
        cout << "====== 语义分析结束 ======\n";
#endif
    }

private:
    ParseTree *tree;
    map<string, SymbolEntry> &symtab;

#ifdef DEBUG
    void logNode(ParseTreeNode *node, const string &stage)
    {
        cout << "\n[语义LOG] 阶段: " << stage << "\n";
        if (!node)
        {
            cout << "  当前结点: <null>\n";
            return;
        }

        cout << "  当前结点: " << node->symbol.name << "\n";
        cout << "  子结点数: " << node->children.size() << "\n";
        for (size_t i = 0; i < node->children.size(); ++i)
        {
            cout << "    ├─ " << node->children[i]->symbol.name << "\n";
        }
    }
#else
    void logNode(ParseTreeNode *, const string &) {}
#endif

    void analyzeNode(ParseTreeNode *node);
    void handleAssignStmt(ParseTreeNode *node);
    void handleIfStmt(ParseTreeNode *node);
    void handleWhileStmt(ParseTreeNode *node);
    void handleDecl(ParseTreeNode *node);

    ExprResult handleArithExpr(ParseTreeNode *node);
    ExprResult handleArithExprPrime(ExprResult inherited, ParseTreeNode *node);
    ExprResult handleMultiExpr(ParseTreeNode *node);
    ExprResult handleMultiExprPrime(ExprResult inherited, ParseTreeNode *node);
    ExprResult handleSimpleExpr(ParseTreeNode *node);

    bool evalBoolExpr(ParseTreeNode *node);

    [[noreturn]] void semanticError(const string &msg);
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
            node->lexeme = tokens[pos].lexeme;
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
InputProcessor::InputProcessor(const string &input)
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

    logNode(node, "analyzeNode");

    const string &name = node->symbol.name;

    if (name == "decl")
        handleDecl(node);
    else if (name == "assgstmt")
        handleAssignStmt(node);
    else if (name == "ifstmt")
        handleIfStmt(node);
    else if (name == "whilestmt")
        handleWhileStmt(node);
    else
    {
        for (auto *child : node->children)
            analyzeNode(child);
    }
}

/* ---------- 语句 ---------- */

void SemanticAnalyzer::handleAssignStmt(ParseTreeNode *node)
{
    logNode(node, "handleAssignStmt");

    ParseTreeNode *idNode = nullptr;
    ParseTreeNode *exprNode = nullptr;

    for (auto *c : node->children)
    {
        if (c->symbol.name == "ID")
            idNode = c;
        else if (c->symbol.name == "arithexpr")
            exprNode = c;
    }

    if (!idNode || !exprNode)
        semanticError("赋值语句结构不完整");

    string var = idNode->lexeme;

    if (!symtab.count(var))
        semanticError("变量未声明: " + var);

    ExprResult rhs = handleArithExpr(exprNode);
    SymbolEntry &entry = symtab[var];

    if (entry.type != rhs.type)
    {
        if (entry.type == ValueType::REAL && rhs.type == ValueType::INT)
        {
            entry.value = rhs.value;
        }
        else
        {
            semanticError("赋值类型不匹配: " + var);
        }
    }
    else
    {
        entry.value = rhs.value;
    }

    entry.hasValue = true;
}

void SemanticAnalyzer::handleIfStmt(ParseTreeNode *node)
{
    logNode(node, "handleIfStmt");

    ParseTreeNode *cond = nullptr;
    ParseTreeNode *thenStmt = nullptr;
    ParseTreeNode *elseStmt = nullptr;

    for (auto *c : node->children)
    {
        if (c->symbol.name == "boolexpr")
            cond = c;
        else if (c->symbol.name == "stmt")
        {
            if (!thenStmt)
                thenStmt = c;
            else
                elseStmt = c;
        }
    }

    if (!cond || !thenStmt)
        semanticError("if 语句结构错误");

    bool condVal = evalBoolExpr(cond);
    if (condVal)
        analyzeNode(thenStmt);
    else if (elseStmt)
        analyzeNode(elseStmt);
}

void SemanticAnalyzer::handleWhileStmt(ParseTreeNode *node)
{
    logNode(node, "handleWhileStmt");

    ParseTreeNode *cond = nullptr;
    ParseTreeNode *body = nullptr;

    for (auto *c : node->children)
    {
        if (c->symbol.name == "boolexpr")
            cond = c;
        else if (c->symbol.name == "stmt")
            body = c;
    }

    if (!cond || !body)
        semanticError("while 语句结构错误");

    while (evalBoolExpr(cond))
        analyzeNode(body);
}

void SemanticAnalyzer::handleDecl(ParseTreeNode *node)
{
    logNode(node, "handleDecl");

    ParseTreeNode *typeNode = nullptr;
    ParseTreeNode *idNode = nullptr;
    ParseTreeNode *valNode = nullptr;

    for (auto *c : node->children)
    {
        if (c->symbol.name == "int" || c->symbol.name == "real")
            typeNode = c;
        else if (c->symbol.name == "ID")
            idNode = c;
        else if (c->symbol.name == "INTNUM" || c->symbol.name == "REALNUM")
            valNode = c;
    }

    if (!typeNode || !idNode)
        semanticError("声明语句结构错误");

    string var = idNode->lexeme;

    if (symtab.count(var))
        semanticError("重复声明变量: " + var);

    ValueType type = (typeNode->symbol.name == "int") ? ValueType::INT : ValueType::REAL;

    SymbolEntry entry;
    entry.type = type;
    entry.hasValue = false;

    if (valNode)
    {
        if (valNode->symbol.name == "INTNUM")
        {
            entry.value = stod(valNode->lexeme);
            entry.hasValue = true;
        }
        else if (valNode->symbol.name == "REALNUM")
        {
            entry.value = stod(valNode->lexeme);
            entry.hasValue = true;
        }
    }

    symtab[var] = entry;
}

/* ---------- 表达式 ---------- */

ExprResult SemanticAnalyzer::handleArithExpr(ParseTreeNode *node)
{
    logNode(node, "handleArithExpr");

    // arithexpr → multexpr arithexprprime
    ExprResult left = handleMultiExpr(node->children[0]);

    // 只有一个 multexpr（无 + -）
    if (node->children.size() == 1)
        return left;

    // arithexprprime
    ParseTreeNode *prime = node->children[1];

    // E
    if (prime->children.size() == 1 &&
        prime->children[0]->symbol.name == "E")
        return left;

    // + multexpr arithexprprime
    string op = prime->children[0]->symbol.name;
    ExprResult right = handleMultiExpr(prime->children[1]);

    ValueType resType = (left.type == ValueType::REAL || right.type == ValueType::REAL) ? ValueType::REAL : ValueType::INT;

    double val = (op == "+") ? left.value + right.value : left.value - right.value;

    ExprResult combined{resType, val};

    // 递归处理后续 arithexprprime
    if (prime->children.size() > 2)
        return handleArithExprPrime(combined, prime->children[2]);

    return combined;
}
ExprResult SemanticAnalyzer::handleArithExprPrime(ExprResult inherited, ParseTreeNode *node)
{
    logNode(node, "handleArithExprPrime");

    // E
    if (node->children.size() == 1 && node->children[0]->symbol.name == "E")
        return inherited;

    string op = node->children[0]->symbol.name;
    ExprResult right = handleMultiExpr(node->children[1]);

    ValueType resType = (inherited.type == ValueType::REAL || right.type == ValueType::REAL) ? ValueType::REAL : ValueType::INT;

    double val = (op == "+") ? inherited.value + right.value : inherited.value - right.value;

    ExprResult combined{resType, val};
    return handleArithExprPrime(combined, node->children[2]);
}
ExprResult SemanticAnalyzer::handleMultiExpr(ParseTreeNode *node)
{
    logNode(node, "handleMultiExpr");

    // multexpr → simpleexpr multexprprime
    ExprResult left = handleSimpleExpr(node->children[0]);

    if (node->children.size() == 1)
        return left;

    return handleMultiExprPrime(left, node->children[1]);
}
ExprResult SemanticAnalyzer::handleMultiExprPrime(ExprResult inherited, ParseTreeNode *node)
{
    logNode(node, "handleMultiExprPrime");

    // E
    if (node->children.size() == 1 && node->children[0]->symbol.name == "E")
        return inherited;

    string op = node->children[0]->symbol.name;
    ExprResult right = handleSimpleExpr(node->children[1]);

    ValueType resType = (inherited.type == ValueType::REAL || right.type == ValueType::REAL) ? ValueType::REAL : ValueType::INT;

    double val = (op == "*") ? inherited.value * right.value : inherited.value / right.value;

    ExprResult combined{resType, val};
    return handleMultiExprPrime(combined, node->children[2]);
}
ExprResult SemanticAnalyzer::handleSimpleExpr(ParseTreeNode *node)
{
    logNode(node, "handleSimpleExpr");

    if (node->symbol.name == "simpleexpr" && node->children.size() == 1)
        return handleSimpleExpr(node->children[0]);

    if (node->symbol.name == "ID")
    {
        string var = node->lexeme;
        if (!symtab.count(var))
            semanticError("使用未声明变量: " + var);

        SymbolEntry &e = symtab[var];
        if (!e.hasValue)
            semanticError("变量未初始化: " + var);

        return {e.type, e.value};
    }

    if (node->symbol.name == "INTNUM")
        return {ValueType::INT, stod(node->lexeme)};

    if (node->symbol.name == "REALNUM")
        return {ValueType::REAL, stod(node->lexeme)};

    semanticError("非法简单表达式: " + node->symbol.name);
}

bool SemanticAnalyzer::evalBoolExpr(ParseTreeNode *node)
{
    logNode(node, "evalBoolExpr");

    ExprResult l = handleArithExpr(node->children[0]);
    ExprResult r = handleArithExpr(node->children[2]);

    ParseTreeNode *opNode = node->children[1];

    if (opNode->children.empty())
        semanticError("布尔运算符结构错误");

    string op = opNode->children[0]->symbol.name;

    if (op == "<")
        return l.value < r.value;
    if (op == "<=")
        return l.value <= r.value;
    if (op == ">")
        return l.value > r.value;
    if (op == ">=")
        return l.value >= r.value;
    if (op == "==")
        return l.value == r.value;
    if (op == "!=")
        return l.value != r.value;

    semanticError("未知布尔运算符: " + op);
}

[[noreturn]] void SemanticAnalyzer::semanticError(const string &msg)
{
    throw logic_error("[语义错误] " + msg);
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

    InputProcessor processor(inputText);
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
#ifdef DEBUG
    cout << "\n====== Program Result ======\n";
#endif
    for (const auto &p : symbolTable)
    {
        cout << p.first << ": " << p.second.value << endl;
    }
}