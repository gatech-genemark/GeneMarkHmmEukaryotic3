//**************************************************
// file: alist.h
//**************************************************

#ifndef LIST_H__
#define LIST_H__

//--------------------------------------------------
template <class T> class ListElement
{
public:
	ListElement( T *t )
	{
		data = t;
		prev = 0;
		next = 0;
	}

	~ListElement()
	{
		delete data;
	}

	ListElement *prev;
	ListElement *next;
	T *data;
};
//--------------------------------------------------
template <class T> class List
{
public:
	List()
	{
		first = 0;
		last = 0;
		current = 0;
	}

	~List()
	{
		Clean();
	}

	virtual void Clean()
	{
		while ( first != 0 ) {
			current = first;
			first = current->next;
			delete current;
		}

		first = 0;
		last = 0;
		current = 0;
	}

	virtual void AddToTail( T *t )
	{
		ListElement<T> *new_e;

		new_e = new ListElement<T>( t );

		if ( first == 0 ) {
			first = new_e;
			last = new_e;
		} else {
			last->next = new_e;
			new_e->prev = last;
			last = new_e;
		}
	}

	virtual void AddToHead( T *t )
	{
		ListElement<T> *new_e;

		new_e = new ListElement<T>( t );

		if ( first == 0 ) {
			first = new_e;
			last = new_e;
		} else {
			first->prev = new_e;
			new_e->next = first;
			first = new_e;
		}
	}

	bool IsEmpty()
	{
		if ( first == NULL )
			return true;
		return false;
	}

	bool IsEnd()
	{
		if ( current == NULL )
			return true;
		return false;
	}

	void SetCurrentToFirst()
	{
		current = first;
	}

	void SetCurrentToLast()
	{
		current = last;
	}

	T *GetCurrent()
	{
		if ( current == 0 )
			return 0;
		return current->data;
	}

	T *GetForward()
	{
		if ( current == 0 )
			return 0;
		else
		{
			ListElement<T> *tmp;

			tmp = current;
			current = current->next;
			return tmp->data;
		}
	}

	T *GetBackward()
	{
		if ( current == 0 )
			return 0;

		ListElement<T> *tmp;

		tmp = current;
		current = current->prev;
		return tmp->data;
	}

protected:

	ListElement<T> *first;
	ListElement<T> *last;
	ListElement<T> *current;
};
//--------------------------------------------------
#endif // LIST_H__

