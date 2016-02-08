/*
 * Copyright (C) 2015 The Queen's University of Belfast
 * Authors:
 *  + Hans Vandierendonck (hvandierendonck@acm.org)
 *
 * Test program for the Swan extension to Cilk (Intel Cilk Plus (TM) version)
 * Based on a test case from github.com/hvdieren/swan
 * 
 * Thie program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * Swan is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with Swan.  If not, see <http://www.gnu.org/licenses/>.
 */

// -*- c++ -*-
#include <cstdlib>
#include <cstring>
#include <cassert>

#include <iostream>

#include <unistd.h>

#include <cilk/cilk.h>
#include <cilk/dataflow.h>

using cilk::versioned;
using cilk::indep;
using cilk::outdep;
using cilk::inoutdep;

void task_in( indep<int> in, int n ) {
    usleep( int(n)*1000 );
    std::cerr << "task_in n=" << int(n) << '\n';
}

void task_out( outdep<int> out, int n ) {
    usleep( int(n)*1000 );
    std::cerr << "task_out n=" << int(n) << '\n';
}

void task_inout( inoutdep<int> inout, int n ) {
    usleep( int(n)*1000 );
    std::cerr << "task_inout n=" << int(n) << '\n';
}

#if OBJECT_COMMUTATIVITY
void task_cinout( cinoutdep<int> inout, int n ) {
    usleep( int(n)*1000 );
    std::cerr << "task_cinout n=" << int(n) << '\n';
}
#endif

void test0( int n ) {
    versioned<int> obj( 0 );

    cilk_spawn  task_inout( (inoutdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i ) {
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    }
    cilk_spawn task_inout( (inoutdep<int>)obj, 0 );

    cilk_sync;
}

void test1( int n ) {
    versioned<int> obj( 0 );

    cilk_spawn task_out( (outdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i ) {
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    }
    cilk_spawn task_out( (outdep<int>)obj, 0 );

    cilk_sync;
}

void test2( int n ) {
    versioned<int> obj( 0 );

    cilk_spawn task_out( (outdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i ) {
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    }
    cilk_spawn task_inout( (inoutdep<int>)obj, 0 );

    cilk_sync;
}

void test3( int n ) {
    versioned<int> obj( 0 );

    cilk_spawn task_out( (outdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i )
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    cilk_spawn task_inout( (inoutdep<int>)obj, 0 );
    for( int i=0; i < n; ++i )
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    cilk_spawn task_inout( (inoutdep<int>)obj, 0 );

    cilk_sync;
}

#if OBJECT_COMMUTATIVITY
void test4( int n ) {
    versioned<int> obj;

    cilk_spawn task_out( (outdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i )
	cilk_spawn task_cinout( (cinoutdep<int>)obj, (n-i)*10 );
    cilk_spawn task_inout( (inoutdep<int>)obj, 0 );
    cilk_spawn task_cinout( (cinoutdep<int>)obj, 0 );
    for( int i=0; i < n; ++i )
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    cilk_spawn task_cinout( (cinoutdep<int>)obj, 0 );

    cilk_sync;
}

void test5( int n ) {
    versioned<int> obj, x;

    cilk_spawn task_out( (outdep<int>)obj, 1000 );
    for( int i=0; i < n; ++i )
	if( i&1 )
	    cilk_spawn task_cinout( (cinoutdep<int>)obj, (n-i)*10 );
	else
	    cilk_spawn task_cinout( (cinoutdep<int>)x, (n-i)*10 );
    for( int i=0; i < n; ++i )
	cilk_spawn task_in( (indep<int>)obj, (n-i)*10 );
    cilk_spawn task_cinout( (cinoutdep<int>)obj, 0 );

    cilk_sync;
}
#endif

int main( int argc, char * argv[] ) {
    if( argc <= 2 ) {
	std::cerr << "Usage: " << argv[0] << " <test> <n>\n";
	return 1;
    }

    int t = atoi( argv[1] );
    int n = atoi( argv[2] );
    switch( t ) {
    case 0:
	test0( n );
	break;
    case 1:
	test1( n );
	break;
    case 2:
	test2( n );
	break;
    case 3:
	test3( n );
	break;
    case 4:
#if OBJECT_COMMUTATIVITY
	test4( n );
#else
	std::cerr << argv[0] << ": commutativity not enabled\n";
#endif
	break;
    case 5:
#if OBJECT_COMMUTATIVITY
	test5( n );
#else
	std::cerr << argv[0] << ": commutativity not enabled\n";
#endif
	break;
    }

    return 0;
}
