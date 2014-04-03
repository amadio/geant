#include <cuda.h>

#ifndef CUDASTD_VECTOR_H
#define CUDASTD_VECTOR_H

namespace cudastd {

  template <typename T> class vector;
  template <typename T> class _iterator;

  template <typename T> class _node {

    friend class vector<T>;
    friend class _iterator<T>;
  
  public:
  private:
    T value;
    _node<T> *next;
    FQUALIFIER _node(const T& t, _node<T> *next) : value(t), next(next) {}
  };

  template <typename T> class _iterator {

    friend class vector<T>;

  public:
    FQUALIFIER T& operator*() { return pointee->value; }
    FQUALIFIER T* operator->() { return &pointee->value; }
    FQUALIFIER const _iterator<T>& operator++() { // pre-increment
      pointee = pointee->next;
      return *this;
    }
    FQUALIFIER const _iterator<T>& operator++(int) { // post-increment
      pointee = pointee->next;
      return *this;
    }
    FQUALIFIER bool operator!=(const _iterator<T>& other) const {
      return this->pointee != other.pointee;
    }

  private:
    _node<T> *pointee;
    FQUALIFIER _iterator(_node<T> *pointee) : pointee(pointee) {}
  };


  template <typename T> class vector {

  public:
    typedef _iterator<T> iterator;

    FQUALIFIER vector() : n(0), head(NULL), tail(NULL) {}
    FQUALIFIER ~vector() { if(head) delete head; }

    FQUALIFIER bool empty() const { return (head == NULL); }
    FQUALIFIER int size() const { return n; }

    FQUALIFIER void clear() {
      if(head != NULL) {
	while(head != tail) {
	  _node<T>* temp = head->next;
	  delete head;
	  head = temp;
	}
	head = NULL;
      }
      n = 0;
    }

    FQUALIFIER T& operator[] (int idx) {
      if(idx >= n || idx < 0) return *end();
      _node<T> *result = head;
      for(int i=0; i<n; i++) {
	if(i==idx) break;
	else result = result->next;
      }
      return result->value;
    }

    FQUALIFIER void push_back(const T& value) {
      _node<T> *newNode = new _node<T>(value, NULL);
      if (head == NULL) {
	head = newNode;
      } else {
	tail->next = newNode;
      }
      tail = newNode;
      n++;
    }

    FQUALIFIER iterator begin() { return _iterator<T>(head); }
    FQUALIFIER iterator end() { return _iterator<T>(NULL); }

  private:
    int n;
    _node<T> *head;
    _node<T> *tail;
  };

}

#endif
