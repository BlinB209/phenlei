inline ActionKeyFUN Task::GetMainAction() { return main_action; };

inline ActionKeyFUN Task::GetPreAction () { return pre_action;  };

inline ActionKeyFUN Task::GetPostAction() { return post_action; };

inline void Task::SetMainAction(ActionKeyFUN action) { main_action = action; };

inline void Task::SetPreAction (ActionKeyFUN action) { pre_action  = action; };

inline void Task::SetPostAction(ActionKeyFUN action) { post_action = action; };

inline ActionKey * Task::GetActionKey() { return actkey; };

inline void Task::SetActionKey(ActionKey *actkey) { this->actkey = actkey; };