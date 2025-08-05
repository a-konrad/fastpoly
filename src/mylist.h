/*------------------------------------------------------------------------*/
/*! \file mylist.h
    \brief contains the class for the doubly linked list.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#include <string>

#ifndef MYLIST_H_
#define MYLIST_H_

class Monom;

class MyList {

	friend class Monom;
	friend class Polynom;

	private: 
	
		// Class used for ListElements.
		class ListElement {
    		public:
        	// Constructor.
        	ListElement(Monom* data)
        	{
        	    this->data = data;
        	    next = NULL;
        	    prev = NULL;
        	}
 
        	// Pointers to next and previous element in the list.
        	ListElement* next;
        	ListElement* prev; 
 
        	// The actual data: Pointers to monomials.
        	Monom* data;
    	};
    	
    // Members of List class.
	ListElement* head;
	ListElement* tail; 
	int size;
	
	public:	
		// Iterator class for the list.
		class Iterator
		{
			protected:
			    // Pointer to current list element.
			    ListElement* iter;
 	
			public:
			    friend MyList;
 		
			    // Constructor.
			    Iterator(void) : iter(NULL) {}
			    Iterator(ListElement* le) : iter(le) {}
			     
			    // Assignment operator.
			    void operator=(ListElement* le) { 
			    	iter = le; 
			    }
 		
			    // Comparison (not equal) operator.
			    bool operator!=(ListElement* le) {
			        if(NULL != iter && iter!=le) return true;
			        else return false;
			    }
		    
			    // Comparison (equal) operator.
			    bool operator==(ListElement* le) {
					if(NULL != iter && iter==le) return true;
					else return false;
    		    }
 	
			    // Incremental operator.
			    void operator++ (int) {
			        if (iter != NULL) iter = iter->next;
			    }
			    
			    void operator-- (int) {
					if(iter != NULL)
					iter = iter->prev;
    		    }
 		
			    // Access to data.
			    Monom* returnData() {   
			        return iter->data;  
			    }
			    
			    ListElement* returnElement() {
			    	return iter;
			    }
		};
	
	public:
		/** Constructor */
		MyList();
		
		/** Destructor */
		virtual ~MyList();
		
		/** Add element to list.

			@param data Monom*

			@return Pointer to new element.
		*/
    	ListElement* add(Monom* data);
    	
    	/** Check whether list is empty.

			@return true if list is empty.
		*/
    	bool isEmpty() { return (head == NULL) ? true : false; }  
    	
    	/** Delete element from list.

			@param data Monom* to delete.
		*/
		void deleteElement(ListElement* element);
		
		/** Delete list. */
		void deleteList();
		
		/** Print list to standard output. */
		void show();
		
		/** Get first element of list.

			@return Pointer to head element.
		*/
		ListElement* begin();
		
		/** Get end of list.

			@return Pointer to end of list.
		*/
		ListElement* end();
		
		/** Get last element of list.

			@return Pointer to tail of list.
		*/
		ListElement* returnTail();	
		
		/** Get size of list.

			@return integer
		*/
		int getSize();

};

#endif /* MYLIST_H_ */
