class control
{
    public:
    int np {};
    int nc {};
    int blockNum {};
    int nodeSetNum {};
    int sideSetNum {};

    double dlt {};
    double timeEnd {};
    double timeCurrent {};
    double timeOutputPeriod {};
    double timetoOutput {};
    int timeStepID {1};
    int timeOutputStepID {1};
    int cellNumMax {};

    
    double winMax {};

    // methods

    control();

    
};