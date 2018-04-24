# author: Ethan Gniot
# This is a simple class for handling the data2 in RunInfo Tables
# from the SRA database. The intention is for users to be able to
# obtain lists of SRA accession numbers by filtering the samples
# based on the metadata in the RunInfo Table. This should be faster
# than filtering the data2 using the GUI provided by the SRA database.

import pandas as pd
import random


class RunTable:
    def __init__(self, run_info_table):
        if isinstance(run_info_table, pd.DataFrame):
            self.df = run_info_table
        else:
            # Read the tab-delimited run_info_table into a pandas Dataframe.
            self.df = pd.read_csv(run_info_table, delimiter='\t')
            # List of the names of the data2 columns.
            self.column_names = list(self.df.columns.values)

    def filter_data(self, column_name, operation, value):
        """
        This method filters out samples in a RunInfo table based on their metadata values.
        :param column_name: string
        :param operation: string
        :param value: scalar
        :return: A RunTable object containing only the samples that meet
        the sorting criteria.
        """
        # Create a dictionary where key=operator as a string, and
        # value=the operator's function. The function will check each value in the
        # column to see if it satisfies the operation. If it does, that value
        # is marked as 'True'. If not, then 'False'. The end result is an array
        # of Trues and Falses that shows which entries in the column satisfy
        # the operation. This array can then be used by other code to
        # select only the samples whose entries were marked as True.

        # Call the operator function from the ops dictionary to check which
        # samples satisfy the sort operation and make a boolean array
        # from the results.
        if isinstance(value, int):
            ops = {
                '<': self.df[column_name].lt(value),
                '>': self.df[column_name].gt(value),
                '==': self.df[column_name].eq(value),
                '<=': self.df[column_name].le(value),
                '>=': self.df[column_name].ge(value),
                '!=': self.df[column_name].ne(value),
            }
            truth_array = ops[operation]
            # Use the boolean array to select only the samples that satisfied
            # the operation (marked "True" in the boolean array).
            filtered_samples = self.df[truth_array]
            # Return the filtered_samples as an instance of the RunTable class
            # so that the user can continue applying filtering methods to the
            # new dataset.
            # Not returning the samples as a RunTable instance would mean the
            # data wouldn't have access to the RunTable class' filtering
            # methods, so the dataset wouldn't be able to be filtered again.
            return RunTable(filtered_samples)

        if isinstance(value, str):
            try:
                ops = {
                    '==': self.df[column_name].eq(value),
                    '!=': self.df[column_name].ne(value)
                }
                truth_array = ops[operation]
                filtered_samples = self.df[truth_array]
                return RunTable(filtered_samples)
            except(TypeError):
                raise TypeError("The data2 type of the metadata may not match "
                                "the data2 type that you entered for the last "
                                "argument ('value').")
            except(KeyError):
                raise KeyError("The 'operation' argument must be either '==' "
                               "or '!=' when the metadata values are not "
                               "numerical.")

    def filter_data_range(self, column_name, in_or_out, value1, value2):
        """
        This method filters out samples in a RunInfo table based on
        ranges of metadata values.

        :param in_or_out: One of two possible strings: "in" or "out". With "in",
         the function will sort data2 and return accession numbers for values
         that fall within the range specified by value1 and value2. "Out" will
         return values that fall outside the range specified by value1
         and value2.
        :param column_name: string
        :param value1: Integer or decimal number.
        :param value2: Integer or decimal number.
        :return A RunTable object containing only the samples that meet
        the sorting criteria.
        """
        try:
            lower = min(value1, value2)
            upper = max(value1, value2)
        except(TypeError, ValueError):
            # TypeError = Raised when an operation or function is applied to an
            # object of inappropriate type. The associated value is a string
            # giving details about the type mismatch.
            # raise AssertionError("Input for 3rd and 4th arguments should be numbers.")
            raise TypeError("Input for 3rd and 4th arguments should be numbers.")

        try:
            assert (in_or_out == 'in') or (in_or_out == 'out')
            if in_or_out == 'in':
                # Get boolean array for samples that fall within the specified
                # bounds of the range.
                truth_array = self.df[column_name].between(lower, upper, inclusive=True)
                # Use boolean array to select the samples.
                filtered_samples = self.df[truth_array]
                # Return accession numbers from the selected samples.
                return RunTable(filtered_samples)
            if in_or_out == 'out':
                col = self.df[column_name]
                # Mask makes cells that evaluate to True NaN.
                # If the function just looks for samples that don't fall inside
                # the specified range, it will also include samples that just
                # have blank entries for that metadata (b/c an entry of "" would
                # not be a number within the specified range, therefore that entry
                # would be marked as 'True"). We only want samples that have
                # metadata entries for the column.
                #

                # To solve this, we can mask the cells that satisfy the
                # specified range and turn all of the values into NaN
                # ("Not a Number")...
                masked_table = col.mask(col.ge(lower) & col.le(upper))  # Needed to use bitwise '&' instead of logical 'and'; still trying to figure out why logical 'and' doesn't work.
                # ... then use the pandas method .isnull() to
                # create a boolean array where we mark entries as True if they
                # are Nan...
                truth_array = masked_table.isnull()
                # ... then use the opposite of the boolean array to select all
                # samples where the entry is NOT an NaN.
                filtered_samples = self.df[-truth_array]
                return RunTable(filtered_samples)
        except(AssertionError):
            raise AssertionError("The 2nd argument (in_or_out) must be either "
                                 "of the strings 'in' or 'out'")

        return None

    def get_accession_numbers(self):
        """
        This method returns an iterable pandas series of sra accession numbers
        for each of the samples in the RunTable.
        :return: A pandas series. Series contains the sra accession numbers for
        each sample.
        """
        acc_nums = self.df['Run']
        return acc_nums

    def random_sample_subset(self, acc_num_list, n):
        """
        This method is for randomly selecting n number of samples from a list
        of accession numbers.
        Function found at "https://stackoverflow.com/questions/2612648/reservoir-sampling".
        :param acc_num_list: list or pandas series; List of Accession nums.
        :param n: int; Number of samples to select from the list.
        :return: list of length n. Each item in the list is a string
        representing the sra accession number of a sequencing run.
        """
        subset = []
        i = 0

        for acc_num in acc_num_list:
            i += 1
            if len(subset) < n:
                subset.append(acc_num)
            else:
                s = int(random.random() * i)
                if s < n:
                    subset[s] = acc_num
        return subset
