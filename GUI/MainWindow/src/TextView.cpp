#include "TextView.h"
#include <QMenu>
#include <QAction>
#include <QContextMenuEvent>
#include <QScrollBar>

#ifdef WIN32
#if defined(VS) && (_MSC_VER >= 1600)
#pragma execution_character_set("utf-8")
#endif
#endif


TextView::TextView(QWidget* parent)
    : QTextEdit(parent)
{
    setReadOnly(true);
	maxCharCount = 100000;
}

void TextView::adjustContents(int maxLength)
{
	QString plainText = toPlainText();
	int _length = plainText.count();
	if (_length > maxLength)
    {
        int position = 0.8*_length;
		if (position > 0)
            plainText.remove(0, position);
        this->setPlainText(plainText);
		moveToEnd();
	}
}

void TextView::moveToEnd()
{
    moveCursor(QTextCursor::End);
    verticalScrollBar()->setValue(verticalScrollBar()->maximum());
}

QString TextView::getText() const
{
    return toPlainText();
}

void TextView::setMaximumCharCount(int max)
{
	maxCharCount = max;
}

void TextView::appendText(const QString &text)
{
    QString temp = toPlainText() + text;
    int _length = temp.count();
    if (_length > maxCharCount)
    {
        int position = 0.8 * _length;
        if (position > 0){
            temp.remove(0, position);
        }
        this->setPlainText(temp);
        moveToEnd();
    }
    else{
        append(text);
    }
}

