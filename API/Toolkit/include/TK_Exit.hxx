template < typename T >
inline void TK_Exit::UnexpectedVarValue(const string &variableName, const T &value)
{
    ostringstream oss;
    oss << "  Error: this situation has not been considered, for " << variableName << " = " << value << endl;
    ExceptionExit(oss.str());
}