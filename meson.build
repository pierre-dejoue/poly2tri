project('poly2tri', 'cpp', default_options : [ 'cpp_std=c++17' ])

#
# Poly2tri lib
#
include = include_directories('poly2tri/include')
lib = static_library('poly2tri', sources : [
	'poly2tri/src/common/node.cpp',
	'poly2tri/src/common/shapes.cpp',
	'poly2tri/src/common/shapes_io.cpp',
	'poly2tri/src/sweep/advancing_front.cpp',
	'poly2tri/src/sweep/cdt.cpp',
	'poly2tri/src/sweep/policy.cpp',
	'poly2tri/src/sweep/sweep.cpp',
	'poly2tri/src/sweep/sweep_context.cpp',
], include_directories : include, cpp_args : '-DP2T_STATIC_EXPORTS')
poly2tri_dep = declare_dependency(include_directories : include, link_with : lib, compile_args : '-DP2T_STATIC_EXPORTS')

#
# Unit tests
#
thread_dep = dependency('threads')
boost_test_dep = dependency('boost', modules : [ 'filesystem', 'unit_test_framework' ], required : false)

if boost_test_dep.found()
	test('Unit Test', executable('unittest_p2t', [
		'unittest/test_cdt.cpp',
		'unittest/test_triangle.cpp',
		'unittest/utility.cpp',
	], dependencies : [poly2tri_dep, boost_test_dep, thread_dep]))
	p2t_all_includes = include_directories('poly2tri/include', 'poly2tri/src')
	test('Unit Test', executable('unittest_predicates', [
		'unittest/test_predicates.cpp'
	], dependencies : [boost_test_dep, thread_dep], include_directories : p2t_all_includes))
endif
