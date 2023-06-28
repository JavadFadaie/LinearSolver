#ifndef AMG_MEASURELIST
#define AMG_MEASURELIST

#include <algorithm>
#include <cassert>
#include <limits>
#include <list>
#include <vector>

template <typename T>
class MeasureList {

    class IndexList;

    struct IndicesOfConstMeasure {

      IndicesOfConstMeasure
        (std::size_t measure
        , IndexList indexList )
      : _measure(measure)
      , _indexList(indexList)
      {};

      std::size_t _measure;
      IndexList   _indexList;
    };

     public:

    MeasureList
      (T  nIndices)
    : _prev(nIndices,0)
    , _next(nIndices,0)
    , _list()
    {}

      void
    insert
      ( T index
      , std::size_t measure ) {

      auto iter = std::find_if
                    (_list.begin()
                    , _list.end()
                    , [&](IndicesOfConstMeasure const & obj)
                      { return measure >= obj._measure; }
                    );

      if(iter == _list.end()) {
        IndexList indexList(_prev,_next);
        indexList.push_back(index);
        _list.emplace_back(measure,indexList);
      }
      else {
        if(iter->_measure == measure) {
          iter->_indexList.push_back(index);
        }
        else {
          IndexList indexList(_prev,_next);
          indexList.push_back(index);
          _list.emplace(iter,measure,indexList);
        }
      }
    }
    
    void
    erase
      ( T index
      , std::size_t measure ) {

      auto iter = std::find_if
                    (_list.begin()
                     , _list.end()
                     , [&](IndicesOfConstMeasure const & obj)
                       { return measure == obj._measure; }
                    );

      if(iter == _list.end()) {
     //   std::stringstream ss;

       std::cout << "Measure " << measure << " does not exist, cannot remove it";
     //  throw std::runtime_error(ss.str());
      }

      iter->_indexList.erase(index);

      if(iter->_indexList.size()==0) {
        _list.erase(iter);
      }

    }

   T
   popIndexWithMaxMeasure() {
      T const retval(_list.front()._indexList.front());
      erase(retval,_list.front()._measure);
      return retval;
    }

    private:

     class IndexList {

      private:

        static T const HEAD = (std::numeric_limits<T>::max()-0);
        static T const TAIL = (std::numeric_limits<T>::max()-1);

      public:

        IndexList
          ( std::vector<T> & prev
          , std::vector<T> & next)
        : _prev(prev)
        , _next(next)
        , _headIndex(HEAD)
        , _tailIndex(TAIL)
        , _size(0)
        {}

        void
        push_back
          (T const index) {
          assert(index != HEAD);
          assert(index != TAIL);

          if(_size == 0) {
            assert(_headIndex == HEAD);
            assert(_tailIndex == TAIL);
            _headIndex = index;
            _tailIndex = index;
            _next[index] = TAIL;
            _prev[index] = HEAD;
            ++_size;
          }
          else {
            _next[_tailIndex] = index;
            _prev[index] = _tailIndex;
            _next[index] = TAIL;
            _tailIndex = index;
            ++_size;
          }
        }
  
       void
        erase(T const index) {

          assert(index != HEAD);
          assert(index != TAIL);
          assert(_size>0);

          if(_size == 1) {
            _headIndex = HEAD;
            _tailIndex = TAIL;
          }
          else
          {
            if(index == _headIndex) {
              _headIndex = _next[index];
              _prev[_next[index]] = HEAD;
            }
            else if(index == _tailIndex) {
              _tailIndex = _prev[index];
              _next[_prev[index]] = TAIL;
            }
            else {
              _next[_prev[index]] = _next[index];
              _prev[_next[index]] = _prev[index];
             }
          }
          --_size;
        }

   
        T
        size() const {
          return _size;
        }

        T
        front() const {
          assert(_size>0);
          return _headIndex;
        }

      private:

        std::vector<T> & _prev;
        std::vector<T> & _next;

        T _headIndex;
        T _tailIndex;
        T _size;

    };

    std::vector<T>                   _prev;
    std::vector<T>                   _next;
    std::list<IndicesOfConstMeasure> _list;
};


#endif


