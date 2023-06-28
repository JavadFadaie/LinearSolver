#ifndef UNIFORMPARTITION_HPP
#define UNIFORMPARTITION_HPP

namespace {

template <typename T>
T div_roundup(T numerator, T denominator)
{ return ( numerator + denominator - static_cast<T>(1) ) / denominator; }

template <class T>
class UniformPartition
{

public:

  UniformPartition( T begin, T end, unsigned short nparts )
  : _begin(begin)
  , _nelems(end-begin)
  , _nparts(nparts)
  {}

  T
  begin( unsigned short ipart) const
  { return _begin
        + div_roundup( static_cast<T>(ipart) * _nelems, _nparts ); }

  T
  end( unsigned short ipart ) const
  { return begin(ipart + 1); }

  T
  nparts() const
  { return _nparts; }

  T
  size(unsigned short ipart ) const
  { return end(ipart) - begin(ipart); }


private:

  T _begin;
  T _nelems;
  T _nparts;

};

}

#endif
