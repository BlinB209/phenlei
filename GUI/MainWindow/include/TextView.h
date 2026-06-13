//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//          PPPPP  H   H  EEEEE  N    N  GGGGG  L      EEEEE  III         +
//          P   P  H   H  E      NN   N  G      L      E       I          +
//          PPPPP  HHHHH  EEEEE  N N  N  G  GG  L      EEEEE   I          +
//          P      H   H  E      N  N N  G   G  L      E       I          +
//          P      H   H  EEEEE  N    N  GGGGG  LLLLL  EEEEE  III         +
//------------------------------------------------------------------------+
//          Platform for Hybrid Engineering Simulation of Flows           +
//          China Aerodynamics Research and Development Center            +
//                     (C) Copyright, Since 2010                          +
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! @file      TextView.h.
//! @brief     TextView 
//! @author    Hleo.
#ifndef TEXTVIEW_H
#define TEXTVIEW_H

#include <QTextEdit>

class TextView : public QTextEdit
{
    Q_OBJECT

public:
    TextView(QWidget* parent = nullptr);

	void adjustContents(int maxLength);
	
	void moveToEnd();

    QString getText()const;

	//Set the maximum number of characters to ensure that too much data is added at one time
    void setMaximumCharCount(int max);

	//Append characters, if the maximum length is exceeded, only the last 20 percent of the text is displayed
    void appendText(const QString& text);

protected:
	
private:
    int maxCharCount;
};

#endif // TEXTVIEW_H
