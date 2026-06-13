template < typename T >
inline void TK_Warning::UnexpectedVarValue(const string &variableName, const T &value)
{
    ostringstream oss;
    oss << "Warning!!! " << "This situation has not been considered, for " << variableName << " = " << value << endl;
    PrintToWindow(oss.str());
    WriteWarningLogFile(oss.str());
}