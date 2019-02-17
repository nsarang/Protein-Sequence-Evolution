namespace utility {


template<class FuncType, class Args>
void Thread_Manager(std::vector<std::function<FuncType> > vecFuncs,
                                    std::vector<Args>& vecDB, bool bVerbose,
                                    std::string sMsg, int nThreads)
{
    std::vector<std::thread> vecThreads;
    int size = vecDB.size(),
        nCount_Now = 0;

    nThreads = (nThreads > 0 ? nThreads : FALLBACK_THRD_NUM);

    for (int i = 0; i < nThreads; ++i) {
        std::vector vecBatch(vecDB.begin() + i * size / nThreads,
                             vecDB.begin() + (i + 1) * (size) / nThreads);
        vecThreads.push_back(std::thread([ =, &nCount_Now] { Processing_Thread(vecFuncs, vecBatch, nCount_Now); })); // BUG #1: must send by copy
    }

    if (bVerbose) {
        utility::Progress_Indicator(sMsg, 0, size);
        while (true) {
            usleep(500 * 1000);
            utility::Progress_Indicator(sMsg, nCount_Now, size);
            if (nCount_Now == size)
                break;
        }
    }

    for (int i = 0; i < nThreads; ++i)
        vecThreads[i].join();
}


template<class FuncType, class Args>
void Processing_Thread(std::vector<std::function<FuncType> > vecFuncs,
                                       std::vector<Args> vecBatch, int& nCount_Now)
{
    static std::mutex _mtx_count;
    
    for (auto& data : vecBatch) {
        for (auto ProcFunction : vecFuncs)
            ProcFunction(data);

        std::lock_guard<std::mutex> lock(_mtx_count);
        nCount_Now++;
    }
}


} // namespace