/*------------------------------------------------------------------------*/
/*! \file mylist.cpp
    \brief contains the class for the doubly linked list.

  Part of FastPoly : A Polynomial Package For Efficient Polynomial Reduction.
  Copyright(C) 2025 Alexander Konrad, University of Freiburg
*/
/*------------------------------------------------------------------------*/

#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include "mylist.h"

MyList::MyList(){
	this->head = this->tail = NULL;
	this->size = 0;
}

MyList::~MyList(){
	this->deleteList();
}

// Add entry.
MyList::ListElement* MyList::add(Monom* data){
	ListElement* newElement = new ListElement(data);
	if(this->isEmpty()) {  // List was empty.
        	this->head = this->tail = newElement;
    } else {  // List not empty. 
      	this->tail->next = newElement;  // Previously last element points to new element.
       	newElement->prev = this->tail;  // New element points to previous last element.
       	this->tail = newElement;  // new Element gets last element.
    }	
  	this->size++;
    return newElement;
}

// Show all list elements.
void MyList::show() {
    std::cout << "Show refList requested. " << std::endl;
    // current node
    ListElement* curr = this->head;
 
    // As long as not reached the end.
    int i = 0;
    while(curr != NULL) {
        // print content
        std::cout << "Element[" << i << "]: " << curr->data << '\n';
 
        // der Nachfolger wird zum aktuellen Knoten
        curr = curr->next;
        i++;
    }
}

void MyList::deleteElement(ListElement* element) {
	if(isEmpty()) return;
	if(element == NULL) {
		std::cout << "Trying to delete not existing element." << std::endl;
		return;
	}
	if(element == this->head) {  // Element to delete is head.
		if (this->head->next==NULL) {  // Only element gets deleted. List is empty afterwards.
			this->head = this->tail = NULL;
		} else {  // Delete head.
			element->next->prev = NULL;
			this->head = element->next;
		}	
	} else if (element == this->tail) {  // Element to delete is tail.
		this->tail = element->prev;
		this->tail->next = NULL; 	
	
	} else {  // Element is somewhere in between.
		element->prev->next = element->next;
		element->next->prev = element->prev;
	}
	delete element;
	this->size--;
}

// Delete complete list and free the memory.
void MyList::deleteList() {
	if(isEmpty()) return;
    	// As long as there are more than one element.
   	while(this->head->next != NULL) {
        	// Take second last element. 
		ListElement* secondLast = this->tail->prev;
        	// Delete last element.
        	delete this->tail;
        	// Make second last element the last element.
        	secondLast->next = NULL;
        	this->tail = secondLast;
    	}
 	// At last delete list head.
 	delete this->head;
 	this->head = NULL;
 	this->size = 0;
}

MyList::ListElement* MyList::begin() {
	return this->head;
}

MyList::ListElement* MyList::end() {
	if (this->tail == NULL) {
		return NULL;
	} else {
		return this->tail->next;
	}
}

MyList::ListElement* MyList::returnTail() {
	return this->tail;
}

int MyList::getSize() {
	return this->size;
}



