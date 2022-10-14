#include "lamda_test.h"
#include "range_test.h"

int 
main(void)
{
	const struct CMUnitTest range_test[] = {
		cmocka_unit_test(valid_range_test),
		cmocka_unit_test(invalid_range_test),
	};

	const struct CMUnitTest lamda_test[] = {
		cmocka_unit_test(node_initialize_and_free_test),
		cmocka_unit_test(node_change_capacity_test),
		cmocka_unit_test(node_append_test),
		cmocka_unit_test(doc_append_test),
		cmocka_unit_test(parse_test),
	};

	int res;

	res = cmocka_run_group_tests(range_test, NULL, NULL);
	if (res != 0)
		return res;

	res = cmocka_run_group_tests(lamda_test, NULL, NULL);
	return res;
}

