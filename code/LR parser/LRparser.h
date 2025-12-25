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

// ========= 实现 =========

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

        Symbol s(tok, SymbolType::TERMINAL);
        if (!grammar.isTerminal(s))
            throw logic_error("Unknown terminal: " + tok);

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

    Symbol s(token, SymbolType::TERMINAL);

    if (!grammar.isTerminal(s))
    {
        throw logic_error("Unknown terminal in input: " + token);
    }

    return s;
}

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

    /* ========= 3. 读取输入 ========= */

    string inputText;
    read_prog(inputText);

    InputProcessor processor(inputText, grammar);
    vector<Token> input = processor.getTokenStream();
}