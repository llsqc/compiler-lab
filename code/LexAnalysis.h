// C语言词法分析器
#include <cstdio>
#include <cstring>
#include <iostream>
#include <map>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
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

/* ========== Token 类 ========== */
struct Token
{
	string lexeme; // 符号内容
	int type;	   // 符号编号

	Token() {}
	Token(const string &s, int t) : lexeme(s), type(t) {}
};

class Lexer
{
private:
	string prog;			 // 源程序
	int pos;				 // 当前扫描位置
	int length;				 // 程序长度
	map<string, int> keyMap; // 符号-编号映射

public:
	// 构造函数
	Lexer(const string &source, const map<string, int> &mp)
		: prog(source), pos(0), length(source.size()), keyMap(mp) {}

	// 是否结束
	bool isEnd() const
	{
		return pos >= length;
	}

	// 获取下一个 Token
	Token nextToken();

	// 词法分析主过程
	void analyze();

private:
	// 跳过空白符
	void skipBlank();

	// 查看当前字符
	char peek() const;

	// 读取当前字符并前进
	char get();

	// 标识符 / 关键字
	Token scanIdentifier();

	// 常数
	Token scanNumber();

	// 字符串
	Token scanString();

	// 注释
	Token scanComment();

	// 运算符 / 界符
	Token scanOperatorOrDelimiter();
};

char Lexer::peek() const
{
	if (isEnd())
		return '\0';
	return prog[pos];
}

char Lexer::get()
{
	if (isEnd())
		return '\0';
	return prog[pos++];
}

void Lexer::skipBlank()
{
	while (!isEnd() && isspace(peek()))
	{
		pos++;
	}
}

Token Lexer::scanIdentifier()
{
	string lex;
	// 第一个字符一定是字母或 _
	lex += get();

	while (!isEnd())
	{
		char ch = peek();
		if (isalnum(ch) || ch == '_')
		{
			lex += get();
		}
		else
		{
			break;
		}
	}

	// 判断是关键字还是标识符
	if (keyMap.count(lex))
	{
		return Token(lex, keyMap[lex]); // 关键字
	}
	return Token(lex, 81); // 标识符
}

Token Lexer::scanNumber()
{
	string lex;
	while (!isEnd() && isdigit(peek()))
	{
		lex += get();
	}
	return Token(lex, 80); // 常数
}

Token Lexer::scanString()
{
	// 当前字符一定是 "
	get();					// 跳过第一个 "
	return Token("\"", 78); // 注意：字符串拆成三部分由外层处理
}

Token Lexer::scanComment()
{
	string lex;
	lex += get(); // '/'
	lex += get(); // '/'

	while (!isEnd() && peek() != '\n')
	{
		lex += get();
	}
	return Token(lex, 79);
}

Token Lexer::scanOperatorOrDelimiter()
{
	// 尝试双字符
	if (pos + 1 < length)
	{
		string two;
		two += prog[pos];
		two += prog[pos + 1];
		if (keyMap.count(two))
		{
			pos += 2;
			return Token(two, keyMap[two]);
		}
	}

	// 单字符
	string one;
	one += get();
	if (keyMap.count(one))
	{
		return Token(one, keyMap[one]);
	}

	return Token(one, -1);
}

Token Lexer::nextToken()
{
	skipBlank();
	if (isEnd())
		return Token("", -1);

	char ch = peek();

	if (isalpha(ch) || ch == '_')
		return scanIdentifier();
	if (isdigit(ch))
		return scanNumber();
	if (ch == '"')
		return scanString();
	if (ch == '/' && pos + 1 < length)
	{
		if (prog[pos + 1] == '/' || prog[pos + 1] == '*')
			return scanComment();
	}

	return scanOperatorOrDelimiter();
}

void Lexer::analyze()
{
	int cnt = 1;
	while (!isEnd())
	{
		skipBlank();
		if (isEnd())
			break;

		Token t = nextToken();
		if (t.type == -1)
			break;

		cout << cnt++ << ": <" << t.lexeme << "," << t.type << ">";
		if (!isEnd())
			cout << endl;
	}
}

void Analysis()
{
	string prog;
	read_prog(prog);
	/* 骚年们 请开始你们的表演 */
	/********* Begin *********/
	map<string, int> keyMap;
	ifstream fin("c_keys.txt");
	string s;
	int id;
	while (fin >> s >> id)
	{
		keyMap[s] = id;
	}
	fin.close();

	Lexer lexer(prog, keyMap);
	lexer.analyze();
	/********* End *********/
}