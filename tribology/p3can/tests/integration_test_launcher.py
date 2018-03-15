import glob
import os
import subprocess
from shutil import copyfile


def backup_current_user_input():
    """Make a copy of the USER_INPUT.py file that can be restored after all integration tests were run"""
    copyfile(os.getcwd() + os.sep + 'USER_INPUT.py', os.getcwd() + os.sep + 'USER_INPUT_BEFORE_TEST.py')


def get_test_input(test_dir):
    """Get a list of all integration test input parameter sets by loading all .py file in the IntegrationTests folder"""
    input_files = sorted(glob.glob('{}{}IntegrationTests{}run{}*.py'.format(test_dir, os.sep, os.sep, os.sep)))
    test_name_list = ('running {} p3can integration tests: {}'.format(len(input_files), len(input_files)*'\n{}'))\
                 .format(*input_files)
    print((max(len(testname) for testname in test_name_list.split('\n'))) * '-')
    print(test_name_list.format(*input_files))
    return input_files


def load_test_input(input_file):
    """Load the current test input"""
    os.chdir('..' + os.sep)
    copyfile(input_file, os.getcwd() + os.sep + 'USER_INPUT.py')
    os.chdir('p3can' + os.sep)


def run_tests(test_inputs):
    """RUn all tests and print error messages to console if something goes wrong"""
    os.chdir('p3can' + os.sep)
    failed_tests = []
    passed_tests = []
    for test_input in test_inputs:
        load_test_input(test_input)
        print('\n{}\nrunning test {}'.format((len(test_input) + 13) * '-', test_input))

        if os.name == 'nt':
            command = [r'python {}{}p3can.py'.format(os.getcwd(), os.sep)]
            cmd = [command]
        else:
            cmd = ['python3 {}/p3can.py'.format(os.getcwd())]
        testrun = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
        out, err = testrun.communicate()

        if err.decode('ascii'):
            err_message = 'p3can integration test failed.'
            print("{}\n{}\n{}".format(len(err_message)*'#', err_message, len(err_message)*'#'))
            try:
                print(out.decode('ascii'), err.decode('ascii'))
            except UnicodeDecodeError:
                print(out.decode('ascii'), err)
            failed_tests.append(test_input)
        else:
            print("p3can integration test passed")
            passed_tests.append(test_input)
        print('{} tests remaining'.format(len(test_inputs)-len(failed_tests)-len(passed_tests)))
    os.chdir('..' + os.sep)
    return passed_tests, failed_tests


def print_summary(passed_tests, failed_tests):
    """Print test summary"""
    print('\n')
    print(max(len(testname) for testname in passed_tests) * '=')
    print(('PASSED TESTS\n{}'.format('{}\n' * len(passed_tests))).format(*passed_tests))
    print(('FAILED TESTS\n{}'.format('{}\n'*len(failed_tests))).format(*failed_tests))


def restore_pre_test_user_input(test_dir):
    """Restore USER_INPUT.txt file to state before running this script"""
    copyfile(os.getcwd() + os.sep + 'USER_INPUT_BEFORE_TEST.py', os.getcwd() + os.sep + 'USER_INPUT.py'.format(os.getcwd()))
    os.remove(os.getcwd() + os.sep + 'USER_INPUT_BEFORE_TEST.py'.format(os.getcwd()))
    os.chdir(test_dir)


def main():
    """Orchestrate test runs"""
    test_dir = os.getcwd()
    os.chdir('..' + os.sep)
    backup_current_user_input()
    test_inputs = get_test_input(test_dir)
    passed_tests, failed_tests = run_tests(test_inputs)
    print_summary(passed_tests, failed_tests)
    restore_pre_test_user_input(test_dir)
    os.chdir(os.getcwd())


if __name__ == "__main__":
    main()
