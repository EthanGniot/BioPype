# author: Ethan Gniot
# This is a simple class for handling the data in RunInfo Tables
# from the SRA database. The intention is for users to be able to
# obtain lists of SRA accession numbers by filtering the samples
# based on the metadata in the RunInfo Table. This should be faster
# than filtering the data using the GUI provided by the SRA database.

import pandas as pd


class RunTable:
    def __init__(self, run_info_table):
        # Read the tab-delimited run_info_table into a pandas Dataframe.
        self.df = pd.read_csv(run_info_table, delimiter='\t')
        # List of the names of the data columns.
        self.column_names = list(self.df.columns.values)

    def sort_data(self, column_name, operation, value):
        """
        This is a function for getting a list of accession numbers for specific
        samples from the RunInfo Table via sorting the table by metadata.
        :param column_name: string
        :param operation: string
        :param value: scalar
        :return: A pandas series of accession numbers for samples that match
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
            true_samples = self.df[truth_array]
            # Return the accession numbers from those selected samples.
            return true_samples['Run']
        if isinstance(value, str):
            try:
                ops = {
                    '==': self.df[column_name].eq(value),
                    '!=': self.df[column_name].ne(value)
                }
                truth_array = ops[operation]
                true_samples = self.df[truth_array]
                return true_samples['Run']
            except(TypeError):
                raise TypeError("The data type of the metadata may not match "
                                "the data type that you entered for the last "
                                "argument ('value').")
            except(KeyError):
                raise KeyError("The 'operation' argument must be either '==' "
                               "or '!=' when the metadata values are not "
                               "numerical.")

    def sort_data_range(self, column_name, in_or_out, value1, value2):
        """
        This is a function for getting a list of accession numbers for specific
        samples from the RunInfo Table based on a range of metadata values.
        :param in_or_out: One of two possible strings: "in" or "out". With "in",
         the function will sort data and return accession numbers for values
         that fall within the range specified by value1 and value2. "Out" will
         return values that fall outside the range specified by value1
         and value2.
        :param column_name: string
        :param value1: Integer or decimal number.
        :param value2: Integer or decimal number.
        """
        try:
            lower = min(value1, value2)
            upper = max(value1, value2)
        except(TypeError, ValueError):
            # TypeError = Raised when an operation or function is applied to an
            # object of inappropriate type. The associated value is a string
            # giving details about the type mismatch.
            raise TypeError("Input for 3rd and 4th arguments should be numbers.")
        try:
            assert (in_or_out == 'in') or (in_or_out == 'out')
            if in_or_out == 'in':
                # Get boolean array for samples that fall within the specified
                # bounds of the range.
                truth_array = self.df[column_name].between(lower, upper, inclusive=True)
                # Use boolean array to select the samples.
                true_samples = self.df[truth_array]
                # Return accession numbers from the selected samples.
                return true_samples['Run']
            if in_or_out == 'out':
                col = self.df[column_name]
                # Mask makes cells that evaluate to True NaN.
                # If the function just looks for samples that don't fall inside
                # the specified range, it will also include samples that just
                # have blank entries for that metadata (b/c an entry of "" would
                # not a number within the specified range, therefore that entry
                # would be marked as 'True"). We only want samples that have
                # metadata entries for the column.
                #

                # To solve this, we can mask the cells that satisfy the
                # specified range and turn all of the values into NaN
                # ("Not a Number")...
                masked_table = col.mask(col.ge(lower) & col.le(upper))
                # ... then use the pandas method .isnull() to
                # create a boolean array where we mark entries as True if they
                # are NaN...
                truth_array = masked_table.isnull()
                # ... then use the opposite of the boolean array to select all
                # samples where the entry is NOT a NaN.
                true_samples = self.df[-truth_array]
                return true_samples['Run']
        except(AssertionError):
            raise AssertionError("The 2nd argument (in_or_out) must be either "
                                 "of the strings 'in' or 'out'")
        return None
